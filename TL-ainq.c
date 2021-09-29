//Standard headers
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/time.h>           //FD_SET, FD_ISSET, FD_ZERO macros

//Crypto headers
#include <time.h>
#include "ecdh_ED25519.h"
#include "randapi.h"


#define TRUE 1
#define FALSE 0
#define PORT 5555

//ED25519 curve parameters
BIG_256_56 q;
ECP_ED25519 P;

//V parameters
char pubvalue[2 * EFS_ED25519 + 1];
octet V = {0, sizeof(pubvalue), pubvalue};

//Definition of a user in AinQ
typedef struct {
    char *ID_I;
    ECP_ED25519 *P_I;
    ECP_ED25519 *R_I;
    BIG_256_56 C_I;
    char *c;
} GL_user;

GL_user GL[500];

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
int gengroupkey(GL_user * GL, csprng *RNG, octet * P_PUB, int cnt)
{
    //Group Key, V and l_k generation
    int pest, i;

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

    for (i = 0; i < cnt; i++){
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

        ECP_ED25519_toOctet(&PI, GL[i].P_I, false);
        ECP_ED25519_toOctet(&RI, GL[i].R_I, false);

        //Hash variables
        char mk[132];
        octet MK = {0, sizeof(mk), mk};

        hash256 htag;
        char hmsg[64];
        octet HMSG = {0, sizeof(hmsg), hmsg};

        //Create the message to be hashed: H_0(ID,R_i,P_i)
        OCT_empty(&MK);
        OCT_jbytes(&MK, GL[i].ID_I, 1);
        OCT_joctet(&MK, &RI);
        OCT_joctet(&MK, &PI);

        //Perform the hashing using SHA256
        SPhash(MC_SHA2, SHA256, &HMSG, &MK);
        //printf("Hash Digest: ");
        //OCT_output(&HMSG);
        //printf("\n");

        //H_0(ID,R_i,P_i).P_pub
        BIG_256_56 hh;
        BIG_256_56_fromBytes(hh, HMSG.val);
        ECP_ED25519_mul(&TP_PUB, hh);
        
        //Point Addition Y_i = R_i + TP_PUB + P_i
        ECP_ED25519_add(&TP_PUB, GL[i].R_I);
        ECP_ED25519_add(&TP_PUB, GL[i].P_I);
        ECP_ED25519_toOctet(&Y_I, &TP_PUB, false); 

        //Compute T_A
        ECP_ED25519_mul(&TP_PUB, l_k);
        ECP_ED25519_toOctet(&T_I, &TP_PUB, false);
        //printf("T_I: ");
        //OCT_output(&T_I);
        //printf("\n"); 

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
        OCT_jbytes(&KK, GL[i].ID_I, 1);
        OCT_joctet(&KK, &RI);
        OCT_joctet(&KK, &PI);

        //Perform the hashing using SHA256
        SPhash(MC_SHA2, SHA256, &KMSG, &KK);
        //printf("XOR Hash Digest: ");
        //OCT_output(&KMSG);
        //printf("\n");

        //Compute XOR to generate Ciphertext
        BIG_256_56 h_val, c_i;
        
        BIG_256_56_fromBytes(h_val, KMSG.val);
        //BIG_256_56_or(c_a, h_val, k_g);   
       
        BIG_256_56_norm(h_val);
        BIG_256_56_norm(k_g);
        for (int j = 0; j < NLEN_256_56; j++)
            c_i[j] = h_val[j] ^ k_g[j];                 //C_I = K_g XOR H_1(Y_A, V, T_A, ID, R_A, P_A)
        
        BIG_256_56_copy(GL[i].C_I, c_i);
        
        GL[i].c = malloc(sizeof(c_i));
        BIG_256_56_toBytes(GL[i].c, c_i);       //Convert Ciphertext to Byte for Sending
        //GL[i].C_I = c_i;
    }

    return 0;
}

//Function display strings
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


//Main Function for implmenting the main script
int main(int argc, char *argv[])
{
	int i, res;
    unsigned long ran;
    char *server_message;

    //Socket Paramters
    int opt = TRUE;
    int master_socket, addrlen, new_socket, client_socket[30], max_clients = 30, activity, valread, sd, max_sd;
    struct sockaddr_in address;

    char buffer[1025];  //data buffer

    fd_set readfds;     //set of socket descriptor

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

	//Generate KGC private and public keys
	gen_secret_value(&RNG, &X, &P_PUB); 

	printf("KGC's Private Key = ");
    OCT_output(&X);
    printf("\n");
    printf("KGC's Public Key = ");
    OCT_output(&P_PUB);
    printf("\n");
    printf("======================================================================================================================\n\n");

	printf("============================================ Begin Socket Creation ====================================================\n");

    //initialize all client sockets to 0
    for (i = 0; i < max_clients; i++){
        client_socket[i] = 0;
    }

    //Create the master socket
    if ((master_socket = socket(AF_INET, SOCK_STREAM, 0)) == 0){
        perror("Master socket creation failed");
        exit(EXIT_FAILURE);
    }

    //Allow master socket to allow multiple connections
    if (setsockopt(master_socket, SOL_SOCKET, SO_REUSEADDR, (char *)&opt, sizeof(opt)) < 0){
        perror("setsockopt");
        exit(EXIT_FAILURE);
    }

    //Socket Address Creation Parameters
    address.sin_family = AF_INET;
    address.sin_addr.s_addr = INADDR_ANY;
    address.sin_port = htons(PORT);

    //bind the socket to the specified port
    if (bind(master_socket, (struct sockaddr *)&address, sizeof(address)) < 0){
        perror("Bind operation failed");
        exit(EXIT_FAILURE);
    }

    //Listen for a maximum of 5 connections
    if (listen(master_socket, 5) < 0){
        perror("listen");
        exit(EXIT_FAILURE);
    }

    //Accept any incoming connections
    addrlen = sizeof(address);
    puts("Waiting for Incoming Connections");

    //count number of connections so far
    int cnt = 0;

    while (TRUE){
        //clear the socket set
        FD_ZERO(&readfds);
        //add the master socket to the set
        FD_SET(master_socket, &readfds);
        max_sd = master_socket;

        //add children sockets to the set
        for (i = 0; i < max_clients; i++){
            //socket descriptor
            sd = client_socket[i];
            //if the descriptor is valid, add to the read list
            if(sd > 0)
                FD_SET(sd, &readfds);
            if(sd > max_sd)
                max_sd = sd;
        }

        //Indefinitely wait for any activity on any of the sockets
        activity = select(max_sd + 1, &readfds, NULL, NULL, NULL);
        if ((activity < 0) && (errno!=EINTR)){
            printf("Error with Selecting Activity\n");
        }

        //Activity of the master socket
        if (FD_ISSET(master_socket, &readfds)){
            //Accept Incoming Connection
            if ((new_socket = accept(master_socket, (struct sockaddr *)&address, (socklen_t *)&addrlen)) < 0){
                perror("Error Accepting Connection");
                exit(EXIT_FAILURE);
            }

            printf("New Connection, Socket FD is %d, IP is: %s, Port: %d\n", new_socket, inet_ntoa(address.sin_addr), ntohs(address.sin_port));

            puts("\n=================== Client Connected =======================");
            //Initial Messages to client
            char *server_message = "WELCOME TEAM MEMBER\n\n========================= Here are the System Public Parameters =====================================\n";
            
            write(new_socket, server_message, strlen(server_message));
            write(new_socket, P_PUB.val, P_PUB.len);
            write(new_socket, p_x, sizeof(p1));
            write(new_socket, p_y, sizeof(p2));
            write(new_socket, q_param, sizeof(q));

            
            //Receive Partial Key Generation Parameters
            char p_i[2 * EFS_ED25519 + 1], r_i[2 * EFS_ED25519 + 1];
            char s_i[2 * EGS_ED25519];
            char clientID[1];
            octet P_I = {0, sizeof(p_i), p_i};
            octet R_I = {0, sizeof(r_i), r_i};
            octet S_I = {0, sizeof(s_i), s_i};
            
            recv(new_socket, clientID, 1, 0);
            printf("\nUser %s's Public Key = ", clientID);

            recv(new_socket, p_i, sizeof(p_i), 0);
            OCT_jbytes(&P_I, p_i, sizeof(p_i));
            OCT_output(&P_I);

            //Generate partial private and public key for the user
            gen_partial_key(&RNG, &P_I, &R_I, &X, &S_I, clientID);
            printf("User %s's Partial Private Key = ", clientID);
            OCT_output(&S_I);
            printf("User %s's Partial Public Key = ", clientID);
            OCT_output(&R_I);
            printf("===============================================================================================================================\n");
            printf("\n");

            //Send the secret parameters back to the user (assumed to be done in a secure manner)
            write(new_socket, S_I.val, S_I.len);
            write(new_socket, R_I.val, R_I.len);

            //Add socket to array of sockets and populate the Group List
            for (i = 0; i < max_clients; i++){
                if(client_socket[i] == 0){
                    client_socket[i] = new_socket;
                    cnt = cnt+1;
                    ECP_ED25519 p_val, r_val;
                    ECP_ED25519_fromOctet(&p_val, &P_I);
                    ECP_ED25519_fromOctet(&r_val, &R_I);

                    //Store Values in the Struct of Arrays for later use
                    GL[i].ID_I = (char *)malloc(strlen(clientID)); //The ID
                    strcpy(GL[i].ID_I, clientID);
                    GL[i].P_I = malloc(sizeof(p_val));  //Partial Public Key from the User
                    ECP_ED25519_copy(GL[i].P_I, &p_val);
                    GL[i].R_I = malloc(sizeof(r_val));  //Partial Public Key from the KGC
                    ECP_ED25519_copy(GL[i].R_I, &r_val);

                    printf("===============================================================================================================================\n\n");
                    break;
                }  
            }

            printf("Updating contents of the Group List\n");
            printf("Number of Users in the Group %d\n\n", cnt);

            //Generating a Group Key for the Group
            printf("\nGenerating Group Key for the Group\n");
            gengroupkey(GL, &RNG, &P_PUB, cnt);
            
            for (i = 0; i < cnt; i++){
                printf("User %s's P_I = ", GL[i].ID_I);
                ECP_ED25519_output(GL[i].P_I);
                printf("User %s's R_I = ", GL[i].ID_I);
                ECP_ED25519_output(GL[i].R_I);
                printf("User %s's C_I = ", GL[i].ID_I);
                BIG_256_56_output(GL[i].C_I);
                printf("\n\n");
            }
            printf("===============================================================================================================================\n");

            //Sending Ciphertext back to the User
            for (i = 0; i < max_clients; i++){
                if(client_socket[i] != 0){
                    printf("Sending User %s's Ciphertext\n", GL[i].ID_I);
                    write(client_socket[i], GL[i].c, 32);
                    write(client_socket[i], V.val, V.len);
                    //break;
                }  
            }

        }

        for (i = 0; i < max_clients; i++){
            sd = client_socket[i];

            if (FD_ISSET( sd , &readfds)){  
                //Check if it was for closing , and also read the 
                //incoming message 
                if ((valread = read( sd , buffer, 1024)) == 0){  
                    //Somebody disconnected , get his details and print 
                    getpeername(sd , (struct sockaddr*)&address , (socklen_t*)&addrlen);  
                    printf("Host disconnected , IP: %s , Port: %d \n" , inet_ntoa(address.sin_addr) , ntohs(address.sin_port));      
                    //Close the socket and mark as 0 in list for reuse 
                    close( sd );  
                    client_socket[i] = 0;  
                }  
                     
                //Echo back the message that came in 
                else{  
                    //set the string terminating NULL byte on the end of the data read 
                    buffer[valread] = '\0';  
                    send(sd , buffer , strlen(buffer) , 0 );  
                }  
            }  
        }
    }
	
    KILL_CSPRNG(&RNG);
	return 0;
}
