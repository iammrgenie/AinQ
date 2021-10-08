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
  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<GpsSubscriber>());
  rclcpp::shutdown();
  return 0;
}
