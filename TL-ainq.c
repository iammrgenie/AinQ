#include <stdio.h>
#include <sys/socket.h>
#include <arpa/inet.h>		//inet_addr
#include <unistd.h>			//write and close
#include <string.h>
#include <pthread.h>		//for multi-threading
#include <stdlib.h>

void *connection_handler(void *);


int main(int argc, char *argv[])
{
	int socket_desc, new_socket, c, *new_sock;
	struct sockaddr_in server, client;
	char *message;

	printf("Begin Socket Creation\n");

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
	server.sin_port = htons(1235);

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
		//Reply to client
		message = "Hello team member, I see you are online now\n";
		write(new_socket, message, strlen(message));

		message = "Assigning you to a handler\n";
		write(new_socket, message, strlen(message));

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
	server_message = "Greetings! This is your handler\n";
	write(sock, server_message, strlen(server_message));

	server_message = "Type something\n";
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