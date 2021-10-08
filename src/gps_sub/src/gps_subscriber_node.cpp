#include <memory>

#include <px4_msgs/msg/sensor_gps.hpp>
#include <chrono>
#include <sstream>
#include <thread>
#include "rclcpp/rclcpp.hpp"

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
     //printf("hello there\n");
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
