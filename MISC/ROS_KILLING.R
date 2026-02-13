k = 0.1

th_ROS_microbe = 0.1
ROS = seq(0,100,1)

ros_kill_factor = ifelse(ROS > th_ROS_microbe,
                         1 / (1 + exp(-k * (ROS - th_ROS_microbe))),
                         0)

plot(ROS, ros_kill_factor)