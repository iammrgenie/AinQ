#include <stdio.h>
#include <sys/socket.h>
#include <arpa/inet.h>		//inet_addr
#include <unistd.h>			//write and close
#include <string.h>
#include <stdlib.h>



int main(int argc, char *argv[])
{
	int socket_desc, new_connection;
	struct sockaddr_in server;
	char client_message[2000], server_reply[2000];

	printf("======================== Begin Socket Creation =========================== \n");

	//create socket
	socket_desc = socket(AF_INET, SOCK_STREAM, 0);

	if (socket_desc == -1)
	{
		printf("Could not create socket\n");
	}

	puts("socket done");

	//initialize connection parameters
	server.sin_addr.s_addr = inet_addr("192.168.50.54");
	server.sin_family = AF_INET;
	server.sin_port = htons(1235);

	//connect to TL server
	if (connect(socket_desc, (struct sockaddr *)&server, sizeof(server)) < 0){
		puts("connection error");
		return 1;
	}

	printf("======================= Connected to Team Leader ======================== \n");

	recv(socket_desc, server_reply, 2000, 0);
	puts(server_reply);
	printf("\n");

	recv(socket_desc, server_reply, 2000, 0);
	puts(server_reply);
	printf("\n");

	for (;;) {
		bzero(client_message, 2000);
		printf("Enter message to the Server \n");
		int n = 0;
		while ((client_message[n++] = getchar()) != '\n')
			;
		write(socket_desc, client_message, sizeof(client_message));
		//receive response from the server
		bzero(server_reply, 2000);
		read(socket_desc, server_reply, sizeof(server_reply));
		printf("Response: %s\n", server_reply);
		if (strncmp("exit", client_message, 4) == 0) {
			printf("Exiting =================== \n");
			break;
		}
	}

	close(socket_desc);
	
	return 0;
}
