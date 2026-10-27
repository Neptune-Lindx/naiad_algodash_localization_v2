from launch import LaunchDescription
from launch.substitutions import PathJoinSubstitution
from launch_ros.actions import Node
from launch_ros.substitutions import FindPackageShare


def generate_launch_description():
    params = PathJoinSubstitution([
        FindPackageShare('naiad_localization_v2'),
        'config',
        'localization_v2.yaml',
    ])
    return LaunchDescription([
        Node(
            package='naiad_localization_v2',
            executable='naiad_localization_v2',
            name='naiad_localization_v2',
            output='screen',
            parameters=[params],
        ),
    ])
