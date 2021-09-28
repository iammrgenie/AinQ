#include <stdio.h>
#include <sys/socket.h>
#include <arpa/inet.h>		//inet_addr
#include <unistd.h>			//write and close
#include <string.h>
#include <stdlib.h>


//Crypto headers
#include <time.h>
#include "ecdh_ED25519.h"
#include "randapi.h"

//ED25519 curve parameters
BIG_256_56 q;
ECP_ED25519 P;

char p_param[2 * EGS_ED25519 +1];
octet P_PARAM = {0, sizeof(p_param), p_param};

//V parameters
char pubvalue[2 * EFS_ED25519 + 1];
octet V = {0, sizeof(pubvalue), pubvalue};

//P_PUB parameters
char ppub[2 * EFS_ED25519 + 1];
octet P_PUB = {0, sizeof(ppub), ppub};


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


int keyretrieval(BIG_256_56 S_i, BIG_256_56 X_i, BIG_256_56 C_i, ECP_ED25519 *P_i, ECP_ED25519 *R_i, char *ID_i, octet *GRPKEY)
{
	//Copy the points P and T into temporary values for use in this function
	ECP_ED25519 V_TP, P_TP;
	ECP_ED25519_copy(&P_TP, &P);
	ECP_ED25519_fromOctet(&V_TP, &V);

	//Y_A and T_A
    char y_i[2 * EFS_ED25519 + 1];
    char t_i[2 * EFS_ED25519 + 1];
    octet Y_i = {0, sizeof(y_i), y_i};
    octet T_i = {0, sizeof(t_i), t_i};

	//Compute T_i = (s_i + x_i).V
	//s_i + x_i
	BIG_256_56 sum, mul_val, v_val;
	BIG_256_56_modadd(sum, S_i, X_i, q);

	//(s_i + x_i) . V
	ECP_ED25519_mul(&V_TP, sum);
	ECP_ED25519_toOctet(&T_i, &V_TP, false); 

    printf("T_I: ");
    OCT_output(&T_i);
    printf("\n");

    //Hash variables
	char kk[65*4 + 1];
	octet KK = {0, sizeof(kk), kk};

	hash256 ktag;
	char kmsg[64];
	octet KMSG = {0, sizeof(kmsg), kmsg};

	//Create the message to be hashed: H_1(Y_i, V, T_i, ID, R_i,P_i)
    char pA1[2 * EFS_ED25519 + 1];
    char pA2[2 * EFS_ED25519 + 1];  
    octet PI = {0, sizeof(pA1), pA1};
    octet RI = {0, sizeof(pA2), pA2};
    ECP_ED25519_toOctet(&PI, P_i, false);
    ECP_ED25519_toOctet(&RI, R_i, false);
    OCT_empty(&KK);
    OCT_joctet(&KK, &V);
    OCT_joctet(&KK, &T_i);
    OCT_jbytes(&KK, ID_i, 1);
    OCT_joctet(&KK, &RI);
    OCT_joctet(&KK, &PI);

    //Perform the hashing using SHA256
	SPhash(MC_SHA2, SHA256, &KMSG, &KK);
	printf("XOR Hash Digest: ");
    OCT_output(&KMSG);
    printf("\n");

    //Compute XOR to generate Ciphertext
    BIG_256_56 h_val, k_val;
    BIG_256_56_fromBytes(h_val, KMSG.val);

 
    BIG_256_56_norm(h_val);
    BIG_256_56_norm(C_i);

    for (int i = 0; i < NLEN_256_56; i++)			// K_g = C_i XOR H_1(V,T,)
        k_val[i] = h_val[i] ^ C_i[i];		
    
    GRPKEY->len = EGS_ED25519;
    BIG_256_56_toBytes(GRPKEY->val, k_val);

    return 0;
}


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
	int socket_desc, new_connection;
	struct sockaddr_in server;
	char client_message[200], server_reply[200];
	char welcome[123];
	char *Q, *p_x, *p_y;
	BIG_256_56 p1, p2;
	Q = malloc(sizeof(q));
	p_x = malloc(sizeof(p1));
	p_y = malloc(sizeof(p2));


	printf("======================== Begin Socket Creation =========================== \n");

	//create socket
	socket_desc = socket(AF_INET, SOCK_STREAM, 0);

	if (socket_desc == -1)
	{
		printf("Could not create socket\n");
	}

	puts("socket done");

	//initialize connection parameters
	server.sin_addr.s_addr = inet_addr("127.0.0.1");
	server.sin_family = AF_INET;
	server.sin_port = htons(5555);

	//connect to TL server
	if (connect(socket_desc, (struct sockaddr *)&server, sizeof(server)) < 0){
		puts("connection error");
		return 1;
	}

	printf("============================================= Connected to Team Leader =========================================================== \n");

	recv(socket_desc, welcome, 123, 0);
	puts(welcome);
	printf("\n");
	bzero(welcome, 123);

	//Retrieve the KGC's Public key and convert to an Octet format
	recv(socket_desc, ppub, sizeof(ppub), 0);
	printf("KGC's Public Key = ");
	OCT_jbytes(&P_PUB, ppub, sizeof(ppub));
	OCT_output(&P_PUB);

	//Retrieve the P parameter and convert to the ED25519 curve type
	recv(socket_desc, p_x, sizeof(q), 0);
	BIG_256_56_fromBytes(p1, p_x);
	recv(socket_desc, p_y, sizeof(q), 0);
	BIG_256_56_fromBytes(p2, p_y);
	ECP_ED25519_set(&P, p1, p2);
	printf("P Parameter = ");
	ECP_ED25519_output(&P);

	//Retrieve the Q parameter and convert to the BIG256 number type
	recv(socket_desc, Q, sizeof(q), 0);
	printf("Q parameter = ");
	BIG_256_56_fromBytes(q, Q);
	BIG_256_56_output(q);
	printf("\n");

	// ======================================================================================================================================================
	//Partial Secret Parameters for entity
	char secA[2 * EGS_ED25519];
    char pubA[2 * EFS_ED25519 + 1];    
    octet X_A = {0, sizeof(secA), secA};
    octet P_A = {0, sizeof(pubA), pubA};

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
	
	/*
    //Receive more miscellaenous messages
	recv(socket_desc, server_reply, 93, 0);
	puts(server_reply);
	bzero(server_reply, 93);
	*/

	//Parameters for Requesting for Partial Key Parameters
	char *ID = "A";
	printf("Key Initialization Process\n");
	

	printf("\nGenerating Partial Secret Parameters\n\n");
	printf("===============================================================================================================================\n");
	//Generate A's private and public keys
	gen_secret_value(&RNG, &X_A, &P_A);

	printf("%s's Partial Private Key = ", ID);
	OCT_output(&X_A);
	printf("\n");
	printf("%s's Partial Public Key = ", ID);
	OCT_output(&P_A);
    printf("\n");
    printf("===============================================================================================================================\n\n");

	//Sending request paramters to the KGC/TL
	write(socket_desc, ID, 1);
	write(socket_desc, P_A.val, P_A.len);

	/*
	for (;;) {
		printf("Welcome to Key Initialization Process\n");
		printf("Enter the identity of the Client.\nEither 'A', 'B', 'C', or 'D'. Enter 'exit' to quit the process\n");
		int n = 0;
		while ((ID[n++] = getchar()) != '\n')
			;

		printf("\nGenerating Partial Secret Parameters\n\n");
		printf("===============================================================================================================================\n");
		//Generate A's private and public keys
		gen_secret_value(&RNG, &X_A, &P_A);

		printf("%c's Partial Private Key = ", ID[1]);
	    OCT_output(&X_A);
	    printf("\n");
	    printf("%c's Partial Public Key = ", ID[1]);
	    OCT_output(&P_A);
	    printf("\n");
	    printf("===============================================================================================================================\n\n");

		//Sending request paramters to the KGC/TL
		write(socket_desc, ID, 1);
		write(socket_desc, P_A.val, P_A.len);

		write(socket_desc, client_message, sizeof(client_message));
		//receive response from the server
		bzero(server_reply, 200);
		read(socket_desc, server_reply, sizeof(server_reply));
		printf("Response: %s\n", server_reply);
		if (strncmp("exit", client_message, 4) == 0) {
			printf("Exiting =================== \n");
			break;
		}
	}
	*/

	close(socket_desc);

	
	return 0;
}
