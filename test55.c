#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/time.h>			//FD_SET, FD_ISSET, FD_ZERO macros

#define TRUE 1
#define FALSE 0
#define PORT 5555

int main (int argc, char *argv[])
{
	int opt = TRUE;
	int master_socket, addrlen, new_socket, client_socket[30], max_clients = 30, activity, i, valread, sd, max_sd;
	struct sockaddr_in address;

	char buffer[1025];  //data buffer

	fd_set readfds;		//set of socket descriptor

	char *message = "Message to be sent";

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

			send(new_socket, message, strlen(message), 0);

			//Add socket to array of sockets
			for (i = 0; i < max_clients; i++){
				if(client_socket[i] == 0){
					client_socket[i] = new_socket;
					printf("Adding Socket to List of Sockets as %d\n", i);
					break;
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

	return 0;
}