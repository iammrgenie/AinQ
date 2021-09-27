//Standard headers
#include <stdio.h>
#include <sys/socket.h>
#include <arpa/inet.h>		//inet_addr
#include <unistd.h>			//write and close
#include <string.h>
#include <pthread.h>		//for multi-threading
#include <stdlib.h>

//Crypto headers
#include <time.h>
#include "ecdh_ED25519.h"
#include "randapi.h"

//ED25519 curve parameters
BIG_256_56 q;
ECP_ED25519 P;

//V parameters
char pubvalue[2 * EFS_ED25519 + 1];
octet V = {0, sizeof(pubvalue), pubvalue};

//Definition of a user in AinQ
typedef struct {
    char *ID_I;
    BIG_256_56 *X_I;
    BIG_256_56 *S_I;
    ECP_ED25519 *P_I;
    ECP_ED25519 *R_I;
    BIG_256_56 *C_I;
    BIG_256_56 *K_I;
} GL_user;

GL_user GL[5000];

//Function to Generate Initial Secret and Public Key pairs
int gen_secret_value(csprng *RNG, octet *SECRET_VALUE, octet *PUBLIC_VALUE) 
{
	int res;
	BIG_256_56 s;
	ECP_ED25519 TMP;

    //Generate private and public key pair using RNG as the seed	
	BIG_256_56_randtrunc(s, q, 2 * CURVE_SECURITY_ED25519, RNG);
    SECRET_VALUE->len = EGS_ED25519;
    BIG_256_56_toBytes(SECRET_VALUE->val, s);

    //Copy the curve generator into a temporary variable
    ECP_ED25519_copy(&TMP, &P);

    //start = clock();

    //EC point multiplication
    ECP_ED25519_mul(&TMP, s);

    //end = clock();
    //cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    //printf("Multiplication took %f seconds to execute \n", cpu_time_used);

	ECP_ED25519_toOctet(PUBLIC_VALUE, &TMP, false);

	res = ECP_ED25519_PUBLIC_KEY_VALIDATE(PUBLIC_VALUE);

	if (res != 0)
    {
        printf("ECP Public Key is invalid!\n");
        return 0;
    }

    return 0;
}

//Function to generate Partial Private and Public Key Pair
int gen_partial_key(csprng *RNG, octet *P_i, octet *R_i, octet *X, octet *PT, char *ID)
{
	int res;
	BIG_256_56 s;
	ECP_ED25519 TMP;

	printf("Partial Key Generation for User %s .... \n\n", ID);

	//r_i
    char secret_value[2 * EGS_ED25519];
    octet SECRET_VALUE = {0, sizeof(secret_value), secret_value};

    //Miscellaenous
	char mk[132];
	octet MK = {0, sizeof(mk), mk};

	hash256 htag;
	char hmsg[64];
	octet HMSG = {0, sizeof(hmsg), hmsg};
    
    //Generate private and public key pair using RNG as the seed	
	BIG_256_56_randtrunc(s, q, 2 * CURVE_SECURITY_ED25519, RNG);		//randomly generate r_i
    SECRET_VALUE.len = EGS_ED25519;
    BIG_256_56_toBytes(SECRET_VALUE.val, s);

    //Copy the curve generator into a temporary variable
    ECP_ED25519_copy(&TMP, &P);

    //EC point multiplication to generate R_i
    ECP_ED25519_mul(&TMP, s);
	ECP_ED25519_toOctet(R_i, &TMP, false); 				//R_i
	res = ECP_ED25519_PUBLIC_KEY_VALIDATE(R_i);
	if (res != 0)
    {
        printf("ECP Public Key is invalid!\n");
        return 0;
    }

    //Create the message to be hashed: H_0(ID,R_i,P_i)
    OCT_empty(&MK);
    OCT_jbytes(&MK, ID, 1);
    OCT_joctet(&MK, R_i);
    OCT_joctet(&MK, P_i);
    printf("Message to be Hashed: \n");
    OCT_output(&MK);
    printf("\n");

    //Perform the hashing using SHA256
	SPhash(MC_SHA2, SHA256, &HMSG, &MK);
	printf("Hash Digest: ");
    OCT_output(&HMSG);
    printf("\n");

    //Modular multiplication xH_0(ID,R_i,P_I) mod q
    BIG_256_56 x_byte, h_byte, mul_byte, sum_byte;
    BIG_256_56_fromBytes(x_byte, X->val);
    BIG_256_56_fromBytes(h_byte, HMSG.val);
    BIG_256_56_modmul(mul_byte, x_byte, h_byte, q);

    //modular addition s_i = r_i + xH_0(ID,R_i,P_I) mod q
    BIG_256_56_modadd(sum_byte, s, mul_byte, q);
    PT->len = EGS_ED25519;
    BIG_256_56_toBytes(PT->val, sum_byte);

    return 0;
}

//Group Key Generation Function
int gengroupkey(GL_user * GL, csprng *RNG, octet * P_PUB)
{
    //Group Key, V and l_k generation
    int pest;

    BIG_256_56 k_g, l_k; //k_g is the group key here
    //Generate random K_g and l_k   
    BIG_256_56_randtrunc(k_g, q, 2 * CURVE_SECURITY_ED25519, RNG);     //random generation of K_g
    BIG_256_56_randtrunc(l_k, q, 2 * CURVE_SECURITY_ED25519, RNG);     //random generation of l_k
    
    ECP_ED25519 TMP;

    //l_k 
    char key[2 * EGS_ED25519];
    octet L = {0, sizeof(key), key};
    L.len = EGS_ED25519;

    printf("Group Key Generated: ");
    BIG_256_56_output(k_g);
    printf("\n\n");
    
    //copy l_k into L
    BIG_256_56_toBytes(L.val, l_k);

    //Compute V
    ECP_ED25519_copy(&TMP, &P);

    ECP_ED25519_mul(&TMP, l_k);
    ECP_ED25519_toOctet(&V, &TMP, false); 
    pest = ECP_ED25519_PUBLIC_KEY_VALIDATE(&V);
    if (pest != 0)
    {
        printf("ECP Public Key is invalid!\n");
        return 0;
    }

    printf("V: ");
    OCT_output(&V);
    printf("\n");

    for (int cnt = 0; cnt < 1500; cnt ++){
        //For each received R_i and P_i, do the following
        //Y_A and T_A points on the curve
        char y_i[2 * EFS_ED25519 + 1];
        char t_i[2 * EFS_ED25519 + 1];
        octet Y_I = {0, sizeof(y_i), y_i};
        octet T_I = {0, sizeof(t_i), t_i};

        ECP_ED25519 TP_PUB;
        ECP_ED25519_fromOctet(&TP_PUB, P_PUB);

        char pA1[2 * EFS_ED25519 + 1];
        char pA2[2 * EFS_ED25519 + 1];  
        octet PI = {0, sizeof(pA1), pA1};
        octet RI = {0, sizeof(pA2), pA2};

        ECP_ED25519_toOctet(&PI, GL[cnt].P_I, false);
        ECP_ED25519_toOctet(&RI, GL[cnt].R_I, false);

        //Hash variables
        char mk[132];
        octet MK = {0, sizeof(mk), mk};

        hash256 htag;
        char hmsg[64];
        octet HMSG = {0, sizeof(hmsg), hmsg};

        //Create the message to be hashed: H_0(ID,R_i,P_i)
        OCT_empty(&MK);
        OCT_jbytes(&MK, GL[cnt].ID_I, 1);
        OCT_joctet(&MK, &RI);
        OCT_joctet(&MK, &PI);

        //Perform the hashing using SHA256
        SPhash(MC_SHA2, SHA256, &HMSG, &MK);
        printf("Hash Digest: ");
        OCT_output(&HMSG);
        printf("\n");

        //H_0(ID,R_i,P_i).P_pub
        BIG_256_56 hh;
        BIG_256_56_fromBytes(hh, HMSG.val);
        ECP_ED25519_mul(&TP_PUB, hh);
        
        //Point Addition Y_i = R_i + TP_PUB + P_i
        ECP_ED25519_add(&TP_PUB, GL[cnt].R_I);
        ECP_ED25519_add(&TP_PUB, GL[cnt].P_I);
        ECP_ED25519_toOctet(&Y_I, &TP_PUB, false); 

        //Compute T_A
        ECP_ED25519_mul(&TP_PUB, l_k);
        ECP_ED25519_toOctet(&T_I, &TP_PUB, false);
        printf("T_I: ");
        OCT_output(&T_I);
        printf("\n"); 

        //Final Hash variables
        char kk[65*4 + 1];
        octet KK = {0, sizeof(kk), kk};

        hash256 ktag;
        char kmsg[64];
        octet KMSG = {0, sizeof(kmsg), kmsg};

        //Create the message to be hashed: H_1(Y_A, V, T_A, ID, R_A,P_A)
        OCT_empty(&KK);
        OCT_joctet(&KK, &V);
        OCT_joctet(&KK, &T_I);
        OCT_jbytes(&KK, GL[cnt].ID_I, 1);
        OCT_joctet(&KK, &RI);
        OCT_joctet(&KK, &PI);

        //Perform the hashing using SHA256
        SPhash(MC_SHA2, SHA256, &KMSG, &KK);
        printf("XOR Hash Digest: ");
        OCT_output(&KMSG);
        printf("\n");

        //Compute XOR to generate Ciphertext
        BIG_256_56 h_val, c_i;
        
        BIG_256_56_fromBytes(h_val, KMSG.val);
        //BIG_256_56_or(c_a, h_val, k_g);   
       
        BIG_256_56_norm(h_val);
        BIG_256_56_norm(k_g);
        for (int i = 0; i < NLEN_256_56; i++)
            c_i[i] = h_val[i] ^ k_g[i];                 //C_I = K_g XOR H_1(Y_A, V, T_A, ID, R_A, P_A)
        
        GL[cnt].C_I = malloc(sizeof(c_i));
        GL[cnt].C_I = &c_i;
    }

    return 0;
}


void *connection_handler(void *);

void displayString(char *sample, int len){
  	int i;
    unsigned char ch;
    for (i = 0; i < len; i++)
    {
        ch = sample[i];
        printf("%02x", ch);
    }
    printf("\n");
}


int main(int argc, char *argv[])
{
	int i, res;
    unsigned long ran;
    char *server_message;

    //socket parameters
    int socket_desc, new_socket, c, *new_sock;
	struct sockaddr_in server, client;
	

	//timer parameters
   	clock_t s1, s2, s3, s4, s5, s6, end;
	double grp_key_time, sec_key_time;

    //Initiate curve parameters from the ROM
    ECP_ED25519_generator(&P);
    BIG_256_56_rcopy(q, CURVE_Order_ED25519);

    printf("\nWelcome to the AinQ Scheme using the ED25519 curve\n\n");
    printf("======================================================================================================================\n\n");
    printf("========================================== System Public Parameters ==================================================\n\n");
    printf("P Parameter = ");
    ECP_ED25519_output(&P);
    printf("\n");
    printf("Q Parameter = ");
    BIG_256_56_output(q);
    printf("\n\n");

    //convert the public parameters to bytes for network transport
    char p_param[2 * EGS_ED25519];
    octet P_PARAM = {0, sizeof(p_param), p_param};
	char *q_param;
	q_param = malloc(sizeof(q));
    BIG_256_56_toBytes(q_param, q);


    BIG_256_56 p1;
    BIG_256_56 p2;
    ECP_ED25519_get(p1, p2, &P);
    char *p_x, *p_y;
    p_x = malloc(sizeof(p1));
    p_y = malloc(sizeof(p2));
    BIG_256_56_toBytes(p_x, p1);
    BIG_256_56_toBytes(p_y, p2);

    //Parameters for KGC
	char seckey[2 * EGS_ED25519];
    char pubkey[2 * EFS_ED25519 + 1];    
    octet X = {0, sizeof(seckey), seckey};
    octet P_PUB = {0, sizeof(pubkey), pubkey};

    //Random number generator parameters
    char raw[100];
    octet RAW = {0, sizeof(raw), raw};
    csprng RNG;                // Crypto Strong RNG

    time((time_t *)&ran);

    RAW.len = 100;              // fake random seed source
    RAW.val[0] = ran;
    RAW.val[1] = ran >> 8;
    RAW.val[2] = ran >> 16;
    RAW.val[3] = ran >> 24;
    for (i = 4; i < 100; i++) RAW.val[i] = i;

    CREATE_CSPRNG(&RNG, &RAW);  // initialise strong RNG

	//start = clock();
	//======================================================================================================================================================

    //======================================================================================================================================================

	//Generate KGC private and public keys
    s3 = clock();
	gen_secret_value(&RNG, &X, &P_PUB);
    s4 = clock();
    sec_key_time = ((double) (s4 - s3)) / CLOCKS_PER_SEC; 

	printf("KGC's Private Key = ");
    OCT_output(&X);
    printf("\n");
    printf("KGC's Public Key = ");
    OCT_output(&P_PUB);
    printf("\n");
     printf("======================================================================================================================\n\n");

	printf("============================================ Begin Socket Creation=====================================================\n");
	//create socket
	socket_desc = socket(AF_INET, SOCK_STREAM, 0);

	if (socket_desc == -1)
	{
		printf("Could not create socket\n");
	}

	puts("socket done");

	//initialize parameters
	server.sin_addr.s_addr = INADDR_ANY;
	server.sin_family = AF_INET;
	server.sin_port = htons(5555);

	//bind address
	if (bind(socket_desc, (struct sockaddr *)&server, sizeof(server)) < 0)
	{
		puts("bind failed");
	}

	puts("bind done");

	//Listen for connections
	listen(socket_desc, 3);

	//Accept connections
	puts("Waiting for incoming connections");
	c = sizeof(struct sockaddr_in);

	while((new_socket = accept(socket_desc, (struct sockaddr *)&client, (socklen_t *)&c)))
	{
		puts("connection accepted");
		//Initial Messages to client
		server_message = "WELCOME TEAM MEMBER\n\n========================= Here are the System Public Parameters =====================================\n";
        printf("Size of welcome message = %lu\n", strlen(server_message));
		write(new_socket, server_message, strlen(server_message));

		write(new_socket, P_PUB.val, P_PUB.len);
		//write(new_socket, P_PARAM.val, P_PARAM.len);
        write(new_socket, p_x, sizeof(p1));
        write(new_socket, p_y, sizeof(p2));
		write(new_socket, q_param, sizeof(q));

		server_message = "===================== Assigning you a Thread ======================\n";
		write(new_socket, server_message, strlen(server_message));

		//initiating threads
		pthread_t sniffer_thread;
		new_sock = malloc(1);
		*new_sock = new_socket;

		if (pthread_create(&sniffer_thread, NULL, connection_handler, (void*) new_sock) < 0){
			perror("could not create thread");
			return 1;
		}

		//Now joing thread
		puts("handler assigned");

	}

	if (new_socket < 0)
	{
		perror("accept failed");
	}

	// =====================================================================================================================================

	close(socket_desc);
	
	return 0;
}

void *connection_handler(void *socket_desc)
{
	//sock descriptor
	int sock = *(int*)socket_desc;
	int read_size;
	char *server_message, client_message[2000];

	//send message to client
	server_message = "================== INITIALIZING COMMUNICATION =======================\n";
	write(sock, server_message, strlen(server_message));

	//receive message from client
	/*
	while ((read_size = recv(sock, client_message, 2000, 0)) > 0){
		//repeat message back to client
		write(sock, client_message, strlen(client_message));
	}

	if(read_size == 0){
		puts("client disconnected");
		fflush(stdout);
	}
	else if(read_size == -1){
		perror("recv failed");
	}
	*/

	for (;;) {
		bzero(client_message, 2000);
		read(sock, client_message, 2000);

		printf("Message from Client: %s\n", client_message);
		server_message = "Well received\n";
		write(sock, server_message, 2000);

		if(strncmp("exit", client_message, 4) == 0) {
			printf("Client Exited\n");
			break;
		}

	}

	//release socket
	free(socket_desc);

	return 0;
}