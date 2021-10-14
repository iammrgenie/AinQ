#include <stdio.h>
#include <sys/socket.h>
#include <sys/un.h>
#include <netinet/in.h>
#include <sys/types.h>
#include <arpa/inet.h>		//inet_addr
#include <unistd.h>			//write and close
#include <string.h>
#include <stdlib.h>
#include <errno.h>

//Crypto headers
#include <time.h>
#include "ecdh_ED25519.h"
#include "randapi.h"

//IPC socket parameters
#define BUF_SIZE 500
#define PORT 6666

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


int keyretrieval(BIG_256_56 S_I, BIG_256_56 X_I, BIG_256_56 C_I, octet *P_I, octet *R_I, char *ID, BIG_256_56 KEY)
{
	//Copy the points P and T into temporary values for use in this function
	ECP_ED25519 V_TP, P_TP;
	ECP_ED25519_copy(&P_TP, &P);
	ECP_ED25519_fromOctet(&V_TP, &V);

	//T_I
    char t_i[2 * EFS_ED25519 + 1];
    octet T_i = {0, sizeof(t_i), t_i};

	//Compute T_i = (s_i + x_i).V
	//s_i + x_i
	BIG_256_56 sum, mul_val, v_val;
	BIG_256_56_modadd(sum, S_I, X_I, q);

	//(s_i + x_i) . V
	ECP_ED25519_mul(&V_TP, sum);
	ECP_ED25519_toOctet(&T_i, &V_TP, false); 

    //printf("T_I: ");
    //OCT_output(&T_i);
    //printf("\n");

    //Hash variables
	char kk[65*4 + 1];
	octet KK = {0, sizeof(kk), kk};

	hash256 ktag;
	char kmsg[64];
	octet KMSG = {0, sizeof(kmsg), kmsg};

	//Create the message to be hashed: H_1(V, T_i, ID, R_i,P_i)
    OCT_empty(&KK);
    OCT_joctet(&KK, &V);
    OCT_joctet(&KK, &T_i);
    OCT_jbytes(&KK, ID, 1);
    OCT_joctet(&KK, R_I);
    OCT_joctet(&KK, P_I);

    //Perform the hashing using SHA256
	SPhash(MC_SHA2, SHA256, &KMSG, &KK);
	//printf("XOR Hash Digest: ");
    //OCT_output(&KMSG);
    //printf("\n");

    //Compute XOR to generate Ciphertext
    BIG_256_56 h_val, k_val;
    BIG_256_56_fromBytes(h_val, KMSG.val);

 
    BIG_256_56_norm(h_val);
    BIG_256_56_norm(C_I);

    for (int j = 0; j < NLEN_256_56; j++)
        k_val[j] = h_val[j] ^ C_I[j];               										// K_g = C_i XOR H_1(V,T,)

    printf("\n======================================================================================================================\n");
    printf("RETREIVED KEY = ");
    BIG_256_56_output(k_val);
    printf("\n======================================================================================================================\n");	

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

void GPS_connect(int TL_address, BIG_256_56 KEY){
	struct sockaddr_in CL_addr;
	int ser_addr, nw_sock, valread;
	int opt = 1;
	int add_len = sizeof(CL_addr);


	if((ser_addr = socket(AF_INET, SOCK_STREAM, 0)) == 0){
		perror("Socket Failed");
		exit(EXIT_FAILURE);
	}

	printf("IPC Socket open at = %d\n", ser_addr);


	if(setsockopt(ser_addr, SOL_SOCKET, SO_REUSEADDR, (char *)&opt, sizeof(opt))){
		perror("Setsocopt");
		exit(EXIT_FAILURE);
	}

	CL_addr.sin_family = AF_INET;
	CL_addr.sin_addr.s_addr = INADDR_ANY;
	CL_addr.sin_port = htons(PORT);

	
	if(bind(ser_addr, (struct sockaddr *) &CL_addr, sizeof(CL_addr)) < 0) {
	    perror("Binding Failed");
	    exit(EXIT_FAILURE);
	}
	
	if(listen(ser_addr, 3) < 0) {
	    perror("Listen Failed");
	    exit(EXIT_FAILURE);
	}

	ssize_t numRead;
	char buf[BUF_SIZE];
	for (;;) {          /* Handle client connections iteratively */
	    printf("Waiting for GPS Data Connection...\n");
	    // NOTE: blocks until a connection request arrives.
	    if((nw_sock = accept(ser_addr, (struct sockaddr *)&CL_addr, (socklen_t *)&add_len)) < 0){
			perror("Accept Failed");
			exit(EXIT_FAILURE);
		} 
	    printf("Connection made at fd = %d\n", nw_sock);

	    // Read at most BUF_SIZE bytes from the socket into buf.
	    while ((numRead = read(nw_sock, buf, BUF_SIZE)) > 0) {
	      // Then, write those bytes from buf into STDOUT.
	      if (write(STDOUT_FILENO, buf, numRead) != numRead) {
	        perror("partial/failed write");
	      }
	    }

	    if (numRead == -1) {
	      perror("read");
	    }

	    if (close(nw_sock) == -1) {
	      perror("close");
    	}
  	}
}


int main(int argc, char *argv[])
{
	int i, res;
    unsigned long ran;
	int socket_desc, new_connection;
	struct sockaddr_in server;
	char client_message[200], server_reply[200];
	char welcome[123];
	BIG_256_56 C_I;
	char *Q, *p_x, *p_y, *c;
	BIG_256_56 p1, p2;
	Q = malloc(sizeof(q));
	p_x = malloc(sizeof(p1));
	p_y = malloc(sizeof(p2));
	c = malloc(sizeof(C_I));


	printf("======================== Begin Socket Creation =========================== \n");

	//create socket
	socket_desc = socket(AF_INET, SOCK_STREAM, 0);

	if (socket_desc == -1)
	{
		printf("Could not create socket\n");
	}

	puts("Socket Created");

	//initialize connection parameters
	server.sin_addr.s_addr = inet_addr("172.18.32.132");
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
	//Security Parameters for the user
	char x_i[2 * EGS_ED25519]; 
	char s_i[2 * EGS_ED25519];
    char p_i[2 * EFS_ED25519 + 1], r_i[2 * EFS_ED25519 + 1];    
    octet X_I = {0, sizeof(x_i), x_i};
    octet P_I = {0, sizeof(p_i), p_i};
    octet S_I = {0, sizeof(s_i), s_i};
    octet R_I = {0, sizeof(r_i), r_i};

    BIG_256_56 x_val, s_val, KEY;

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
	

	//Parameters for Requesting for Partial Key Parameters
	char ID;
	
	for (;;) {
		printf("\n\nWelcome to Key Initialization Process\n");
		printf("Enter the identity of the Client.\nEither 'A', 'B', 'C', or 'D'. Enter 'exit' to quit the process\n");
		scanf("%c", &ID);

		printf("\nGenerating Partial Secret Parameters for User %c ........\n\n", ID);
		printf("===============================================================================================================================\n");
		//Generate A's private and public keys
		gen_secret_value(&RNG, &X_I, &P_I);

		printf("%c's Secret Key = ", ID);
	    OCT_output(&X_I);
	    BIG_256_56_fromBytes(x_val, X_I.val);
	    printf("\n");
	    printf("%c's Partial Public Key = ", ID);
	    OCT_output(&P_I);
	    printf("\n");
	    printf("===============================================================================================================================\n\n");

		//Sending request paramters to the KGC/TL
		write(socket_desc, &ID, 1);
		write(socket_desc, P_I.val, P_I.len);

		//Receive partial parameters from the KGC/TL
		recv(socket_desc, s_i, 32, 0);	//S_I value
		BIG_256_56_fromBytes(s_val, s_i);
		printf("Partial Private Key from KGC = ");
		OCT_jbytes(&S_I, s_i, 32);
		OCT_output(&S_I);

		recv(socket_desc, r_i, sizeof(r_i), 0);	//R_I value
		printf("Partial Public Key from KGC = ");
		OCT_jbytes(&R_I, r_i, sizeof(r_i));
		OCT_output(&R_I);

		//Retrieve the Ciphertext and V value from the TL and convert to the BIG256 number type
		recv(socket_desc, c, 32, 0);
		printf("Received Ciphertext = ");
		BIG_256_56_fromBytes(C_I, c);
		BIG_256_56_output(C_I);
		printf("\n");
		recv(socket_desc, pubvalue, sizeof(pubvalue), 0);
		printf("V value from KGC = ");
		OCT_jbytes(&V, pubvalue, sizeof(pubvalue));
		OCT_output(&V);

		//Retrieve the Generated Key
		keyretrieval(s_val, x_val, C_I, &P_I, &R_I, &ID, KEY);
		//printf("Retrieved Key = ");
		//BIG_256_56_output(KEY);

		GPS_connect(socket_desc, KEY);

		/*
		recv(socket_desc, c, sizeof(C_I), 0);
		printf("Updated Ciphertext = ");
		BIG_256_56_fromBytes(C_I, c);
		BIG_256_56_output(C_I);
		printf("\n");
		recv(socket_desc, pubvalue, sizeof(pubvalue), 0);
		printf("Updated V value from KGC = ");
		OCT_jbytes(&V, pubvalue, sizeof(pubvalue));
		OCT_output(&V);
		keyretrieval(s_val, x_val, C_I, &P_I, &R_I, &ID, KEY);

		recv(socket_desc, c, sizeof(C_I), 0);
		printf("Updated Ciphertext = ");
		BIG_256_56_fromBytes(C_I, c);
		BIG_256_56_output(C_I);
		printf("\n");
		recv(socket_desc, pubvalue, sizeof(pubvalue), 0);
		printf("Updated V value from KGC = ");
		OCT_jbytes(&V, pubvalue, sizeof(pubvalue));
		OCT_output(&V);
		keyretrieval(s_val, x_val, C_I, &P_I, &R_I, &ID, KEY);
		*/

		//receive response from the server
		bzero(server_reply, 200);
		read(socket_desc, server_reply, sizeof(server_reply));
		printf("Response: %s\n", server_reply);
		if (strncmp("exit", client_message, 4) == 0) {
			printf("Exiting =================== \n");
			break;
		}
	}

	close(socket_desc);

	KILL_CSPRNG(&RNG);
	return 0;
}
