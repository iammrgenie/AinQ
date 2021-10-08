#include <memory>

#include <sys/types.h>
#include <sys/socket.h>
#include <sys/un.h>

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

using namespace std::chrono;
using std::placeholders::_1;

class GpsSubscriber : public rclcpp::Node
{
  public:
    GpsSubscriber()
    : Node("gps_subscriber")
    {      
      subscription_ = this->create_subscription<px4_msgs::msg::SensorGps>(
      "/sad01/SensorGps_PubSubTopic", 10, std::bind(&GpsSubscriber::topic_callback, this, _1));
      printf("Hello Hi\n");
    }

  private:
    void topic_callback(const px4_msgs::msg::SensorGps::SharedPtr msg) const
    {
     RCLCPP_INFO(this->get_logger(), "Coordinates = Latitude: %d Longitude: %d Altitude: '%d'", msg->lat, msg->lon, msg->alt);
    }
    rclcpp::Subscription<px4_msgs::msg::SensorGps>::SharedPtr subscription_;
};

int main(int argc, char * argv[])
{
  struct sockaddr_un addr;
  ssize_t numRead;
  char buf[BUF_SIZE];

  // Create the node socket with domain: AF_UNIX, type: SOCK_STREAM, protocol: 0
  int sfd = socket(AF_UNIX, SOCK_STREAM, 0);
  printf("Client socket fd = %d\n", sfd);

  if (sfd == -1) {
    perror("socket");
  }

  memset(&addr, 0, sizeof(struct sockaddr_un));
  addr.sun_family = AF_UNIX;
  strncpy(addr.sun_path, SV_SOCK_PATH, sizeof(addr.sun_path) - 1);

  if (connect(sfd, (struct sockaddr *) &addr, sizeof(struct sockaddr_un)) == -1) {
    perror("connect");
  }

  // Read at most BUF_SIZE bytes from STDIN into buf.
  //while ((numRead = read(STDIN_FILENO, buf, BUF_SIZE)) > 0) {
    // Then, write those bytes from buf into the socket.
    //if (write(sfd, buf, numRead) != numRead) {
      //perror("partial/failed write");
    //}
  //}

  //if (numRead == -1) {
    //perror("read");
  //}

  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<GpsSubscriber>());
  rclcpp::shutdown();

  // Closes our socket; server sees EOF.
  exit(EXIT_SUCCESS);
  return 0;
}
