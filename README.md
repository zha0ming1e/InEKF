# InEKF 

**InEKF** is a light-weight library for **In**variant **E**xtended **K**alman **F**ilter on matrix Lie groups. InEKF is implemented for systems whose states evolve on the abstract $\mathrm{N}$-dimensional $\mathrm{K}$ direct isometry matrix Lie groups $\mathrm{SE_K(N)}$ (including the 3-dimensional special orthogonal group $\mathrm{SO(3)}$, the 3-dimensional $\mathrm{K}$ direct isometry special Euclidean group $\mathrm{SE_K(3)}$ and $\mathrm{N}$-dimensional translation matrix group $\mathrm{T(N)}$) and Euclidean vector spaces. Except for InEKF, the library can be extended to other Kalman filters, e.g. the conventional Extended Kalman Filter (EKF) or the Iterated Error-State Extended Kalman Filter (IESEKF) for multiple robot platforms, such as legged robots and Unmanned Aerial Vehicles (UAVs). 



## Build 

### Prerequisites 

- **GCC** and [**CMake**](https://cmake.org/) 

- [**Eigen3**](http://eigen.tuxfamily.org/): Linear algebra 

    On Ubuntu 18.04, we can install CMake and Eigen3 following: 

    ```bash
    # gcc and cmake 
    sudo apt-get install gcc cmake 
    # Eigen3 
    sudo apt-get install libeigen3-dev
    ```

### Build the Library 

- We can build the InEKF library following: 

    ```bash
    git clone https://github.com/zha0ming1e/InEKF 
    cd ./InEKF/ 
    mkdir build && cd build/ 
    cmake .. && make -j8 
    # or: cmake .. && make -j$(( $( cat /proc/cpuinfo | grep "processor" | sort -u | wc -l ) - 1 ))
    ```

### Include InEKF in Your Projects 

- InEKF can be included in your cmake project by adding the InEKF directory to your CMakeLists.txt: 

    ```bash
    find_package(inekf) 
    include_directories(${inekf_INCLUDE_DIRS}) 
    ```



## References 

- Barrau, A. and Bonnabel, S., 2016. The invariant extended Kalman filter as a stable observer. *IEEE Transactions on Automatic Control*, *62*(4), pp.1797-1812. 

    ```
    @article{barrau2016invariant,
      title={The invariant extended Kalman filter as a stable observer},
      author={Barrau, Axel and Bonnabel, Silv{\`e}re},
      journal={IEEE Transactions on Automatic Control},
      volume={62},
      number={4},
      pages={1797--1812},
      year={2016},
      publisher={IEEE}
    }
    ```

- Hartley, R., Ghaffari, M., Eustice, R.M. and Grizzle, J.W., 2020. Contact-aided invariant extended Kalman filtering for robot state estimation. *The International Journal of Robotics Research*, *39*(4), pp.402-430. 

    ```
    @article{hartley2020contact,
      title={Contact-aided invariant extended Kalman filtering for robot state estimation},
      author={Hartley, Ross and Ghaffari, Maani and Eustice, Ryan M and Grizzle, Jessy W},
      journal={The International Journal of Robotics Research},
      volume={39},
      number={4},
      pages={402--430},
      year={2020},
      publisher={SAGE Publications Sage UK: London, England}
    }
    ```

- Bloesch, M., Hutter, M., Hoepflinger, M.A., Leutenegger, S., Gehring, C., Remy, C.D. and Siegwart, R., 2013. State estimation for legged robots-consistent fusion of leg kinematics and IMU. *Robotics*, *17*, pp.17-24. 

    ```
    @article{bloesch2013state,
      title={State estimation for legged robots-consistent fusion of leg kinematics and IMU},
      author={Bloesch, Michael and Hutter, Marco and Hoepflinger, Mark A and Leutenegger, Stefan and Gehring, Christian and Remy, C David and Siegwart, Roland},
      journal={Robotics},
      volume={17},
      pages={17--24},
      year={2013},
      publisher={MIT Press}
    }
    ```

- Ramadoss, P., Romualdi, G., Dafarra, S., Chavez, F.J.A., Traversaro, S. and Pucci, D., 2021, May. Diligent-kio: A proprioceptive base estimator for humanoid robots using extended kalman filtering on matrix lie groups. In *2021 IEEE International Conference on Robotics and Automation (ICRA)* (pp. 2904-2910). IEEE. 

    ```
    @inproceedings{ramadoss2021diligent,
      title={Diligent-kio: A proprioceptive base estimator for humanoid robots using extended kalman filtering on matrix lie groups},
      author={Ramadoss, Prashanth and Romualdi, Giulio and Dafarra, Stefano and Chavez, Francisco Javier Andrade and Traversaro, Silvio and Pucci, Daniele},
      booktitle={2021 IEEE International Conference on Robotics and Automation (ICRA)},
      pages={2904--2910},
      year={2021},
      organization={IEEE}
    }
    ```

- Ramadoss, P., Romualdi, G., Dafarra, S., Traversaro, S. and Pucci, D., 2022, October. Comparison of EKF-Based Floating Base Estimators for Humanoid Robots with Flat Feet. In *2022 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS)* (pp. 6780-6787). IEEE. 

    ```
    @inproceedings{ramadoss2022comparison,
      title={Comparison of EKF-Based Floating Base Estimators for Humanoid Robots with Flat Feet},
      author={Ramadoss, Prashanth and Romualdi, Giulio and Dafarra, Stefano and Traversaro, Silvio and Pucci, Daniele},
      booktitle={2022 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS)},
      pages={6780--6787},
      year={2022},
      organization={IEEE}
    }
    ```