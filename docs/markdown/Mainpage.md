# Sandals

`Astro` is a library for orbital mechanics and spacecraft dynamics. This is a refactored version of the [original library](https://github.com/ebertolazzi/Astro) developed by Enrico Bertolazzi (University of Trento, Department of Industrial Engineering).

## Installation

### Prerequisites

`Astro` carries a set of minimal dependencies. Here's what you need to get started:

- C++ compiler with [`C++17`](https://en.cppreference.com/w/cpp/17) support
- [`CMake`](https://cmake.org) >= 3.10
- [`Eigen3`](https://eigen.tuxfamily.org/index.php?title=Main_Page) >= 3.4.0
- [`Matplot++`](https://alandefreitas.github.io/matplotplusplus/) >= 1.2.0 for plotting (optional).
- [`Sandals`](https://StoccoDavide.github.io/Sandals) ODEs/DAEs C++17 integrator library (optional)

The `Matplot++` library is optional and only required if you want to plot the results of your simulations and is currently used to plot the results of some tests. Similarly, the `Sandals` library is optional and only required if you want to integrate ODEs/DAEs in your simulations.

### Build

1. Clone the repository.
   ```bash
   git clone git@github.com:StoccoDavide/Astro.git
    ```
2. Create a build directory.
    ```bash
    mkdir build
    cd build
    ```
3. Build the project.
   ```bash
    cmake ..
    ```
4. Install the project.
    ```bash
    make install
    ```

## Authors

- Davide Stocco <br>
  University of Trento <br>
  Department of Industrial Engineering <br>
  email: davide.stocco@unitn.it

- Enrico Bertolazzi <br>
  University of Trento <br>
  Department of Industrial Engineering <br>
  email: enrico.bertolazzi@unitn.it

Aka...

```
▗▄▄▄  ▄   ▄  ▐▌    ▗▞▀▜▌▄▄▄▄     ▐▌    ▗▄▄▖ ▗▞▀▚▖ ▄▄▄ ▄   ▄
▐▌  █ █   █  ▐▌    ▝▚▄▟▌█   █    ▐▌    ▐▌ ▐▌▐▛▀▀▘█    █   █
▐▌  █  ▀▄▀▗▞▀▜▌         █   █ ▗▞▀▜▌    ▐▛▀▚▖▝▚▄▄▖█     ▀▀▀█
▐▙▄▄▀     ▝▚▄▟▌               ▝▚▄▟▌    ▐▙▄▞▘          ▄   █
                                                       ▀▀▀
```

## License

The `Astro` project is distributed under the BSD 2-Clause License - see the [LICENSE](https://StoccoDavide.github.io/Astro/LICENSE) file for details.

Here's what the license entails:

1. Anyone can copy, modify and distribute this software.
2. You have to include the license and copyright notice with each and every distribution.
3. You can use this software privately.
4. You can use this software for commercial purposes.
5. This software is provided without warranty.
6. The software author or license can not be held liable for any damages inflicted by the software.
