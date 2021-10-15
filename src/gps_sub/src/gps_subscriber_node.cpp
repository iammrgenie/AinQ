#include <memory>

#include <sys/types.h>
#include <sys/socket.h>
#include <sys/un.h>
#include <netinet/in.h>
#include <sys/types.h>
#include <arpa/inet.h>    //inet_addr
#include <errno.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <px4_msgs/msg/sensor_gps.hpp>
#include <chrono>
#include <sstream>
#include <thread>
#include "rclcpp/rclcpp.hpp"

#define BUF_SIZE 500
#define SV_SOCK_PATH "tpf_unic_sock.server"
#define PORT 6666

using namespace std::chrono;
using std::placeholders::_1;

int node_sock = 0;

class GpsSubscriber : public rclcpp::Node
{
  public:
    GpsSubscriber()
    : Node("gps_subscriber")
    {      
      subscription_ = this->create_subscription<px4_msgs::msg::SensorGps>(
      "/sad01/SensorGps_PubSubTopic", 10, std::bind(&GpsSubscriber::topic_callback, this, _1));
    }

  private:
    void topic_callback(const px4_msgs::msg::SensorGps::SharedPtr msg) const
    {
     RCLCPP_INFO(this->get_logger(), "Latitude Value = %d ", msg->lat);
     char *lat_val;
     std::sprintf(lat_val, "%d", msg->lat);
     printf("Size of Latitude = %ld \n", strlen(lat_val));
     //write(node_sock, lat_val, strlen(lat_val));
    }
    rclcpp::Subscription<px4_msgs::msg::SensorGps>::SharedPtr subscription_;
};

int main(int argc, char * argv[])
{
  struct sockaddr_in N_addr;

  ssize_t numRead;
  //char buf[BUF_SIZE];

  // Create a new client socket with domain: AF_UNIX, type: SOCK_STREAM, protocol: 0
  if((node_sock = socket(AF_INET, SOCK_STREAM, 0)) < 0){
    perror("Socket Failed");
    exit(EXIT_FAILURE);
  }

  printf("Node socket at = %d\n", node_sock);

  N_addr.sin_addr.s_addr = inet_addr("127.0.0.1");
  N_addr.sin_family = AF_INET;
  N_addr.sin_port = htons(PORT);

  /*
  // Connects the active socket via the listening socket
  if(connect(node_sock, (struct sockaddr *)&N_addr, sizeof(N_addr)) < 0){
    perror("Connection Error");
    exit(EXIT_FAILURE);
  }
  */

  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<GpsSubscriber>());
  rclcpp::shutdown();

  // Closes our socket; server sees EOF.
  exit(EXIT_SUCCESS);
  return 0;
}
