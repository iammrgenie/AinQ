#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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


int main()
{
    int i, res;
    unsigned long ran;

   	clock_t s1, s2, s3, s4, s5, s6, end;
	double grp_key_time, sec_key_time;

    //Initiate curve parameters from the ROM
    ECP_ED25519_generator(&P);
    BIG_256_56_rcopy(q, CURVE_Order_ED25519);

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

	printf("\nWelcome to the AinQ Scheme using the ED25519 curve\n\n");
	printf("===============================================================================================================================\n");

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
    printf("===============================================================================================================================\n");

    //Iteration for Performance purposes
    //Secret Key, Public key, and Partial Key Generation Processes
    for (int i = 0; i < 1500; i ++){
		char ID_I[5];
		int tt = snprintf(ID_I, 5, "%d", i);
    	char pubA1[2 * EFS_ED25519 + 1];
    	char pubA2[2 * EFS_ED25519 + 1];  
    	char xA[2 * EGS_ED25519];
    	char pA[2 * EGS_ED25519];
    	octet S_I = {0, sizeof(pA), pA};   
    	octet X_I = {0, sizeof(xA), xA};
    	octet P_I = {0, sizeof(pubA1), pubA1};
    	octet R_I = {0, sizeof(pubA2), pubA2};

    	//Generate private and public key pair for A
    	gen_secret_value(&RNG, &X_I, &P_I);
    	printf("User %d's Secret Value = ", i);
    	OCT_output(&X_I);
    	printf("\n");
    	printf("User %d's Public Value = ", i);
    	OCT_output(&P_I);
    	printf("\n");
    	printf("===============================================================================================================================\n");
    	printf("\n");

    	//Generate partial private and public key for A
		gen_partial_key(&RNG, &P_I, &R_I, &X, &S_I, ID_I);
		printf("User %d's Partial Private Key = ", i);
		OCT_output(&S_I);
		printf("User %d's Partial Public Key = ", i);
		OCT_output(&R_I);
    	printf("===============================================================================================================================\n");
    	printf("\n");

    
    	BIG_256_56 s_val, x_val;
    	ECP_ED25519 p_val, r_val;

    	BIG_256_56_fromBytes(x_val, X_I.val);
    	BIG_256_56_fromBytes(s_val, S_I.val);
    	ECP_ED25519_fromOctet(&p_val, &P_I);
    	ECP_ED25519_fromOctet(&r_val, &R_I);

    	//Store Values in the Struct of Arrays for later use
    	GL[i].ID_I = (char *)malloc(strlen(ID_I));
    	strcpy(GL[i].ID_I, ID_I);

    	GL[i].X_I = malloc(sizeof(x_val));
    	GL[i].X_I = &x_val;

    	GL[i].S_I = malloc(sizeof(s_val));
    	GL[i].S_I = &s_val;

    	GL[i].P_I = malloc(sizeof(p_val));
    	ECP_ED25519_copy(GL[i].P_I, &p_val);

    	GL[i].R_I = malloc(sizeof(r_val));
    	ECP_ED25519_copy(GL[i].R_I, &r_val);

    	printf("===============================================================================================================================\n");
    	printf("\n");

        /*
    	char c_i[2 * EGS_ED25519];
    	octet C_I = {0, sizeof(c_i), c_i};
    	char k_i[2 * EGS_ED25519];
		octet K_I = {0, sizeof(k_i), k_i};
        */
	}

    s1 = clock();
    //User Ciphertext Generation
    gengroupkey(GL, &RNG, &P_PUB);
    s2 = clock();
    grp_key_time = ((double) (s2 - s1)) / CLOCKS_PER_SEC; 
    

    /*
    printf("===============================================================================================================================\n");
    printf("\n");

    for (int j = 0; j < 1500; j++){
        //Group Key Retrieval
        char gkey[2 * EGS_ED25519];
        octet GRPKEY = {0, sizeof(gkey), gkey};

        printf("C_I: ");
        BIG_256_56_output(GL[j].C_I);
        printf("\n");

        keyretrieval(GL[j].S_I, GL[j].X_I, GL[j].C_I, GL[j].P_I, GL[j].R_I, GL[j].ID_I, &GRPKEY);
        printf("User %s's Retrieved Group Key's Value = ", GL[j].ID_I);
        OCT_output(&GRPKEY);
        printf("\n");
        printf("===============================================================================================================================\n");
        printf("\n");

        BIG_256_56 k_val;
        BIG_256_56_fromBytes(k_val, GRPKEY.val);
        GL[i].K_I = malloc(sizeof(k_val));
        GL[i].K_I = &k_val;
        //printf("Done\n");
    }

    */

    printf("Total execution time for the Secret Key Generation = %f seconds\n", sec_key_time);
    printf("Total execution time for Group Key Generation = %f seconds\n", grp_key_time);

    /*
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("Total took %f seconds to execute \n", cpu_time_used);
    */
    /*
    for (int z = 0; z < 1; z ++)
    {
    	printf("User %s's Group Key = ", GL[z].ID_I);
    	BIG_256_56_output(GL[z].K_I);
    	printf("\n");
    }
    */

    

    KILL_CSPRNG(&RNG);
    return 0;
}