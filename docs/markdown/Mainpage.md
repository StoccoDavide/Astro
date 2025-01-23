# Sandals

`Astro` is a library for orbital mechanics and spacecraft dynamics. This is a refactored version of the [original library](https://github.com/ebertolazzi/Astro) developed by Enrico Bertolazzi (University of Trento, Department of Industrial Engineering).

## Installation

### Quick and dirty

`Astro` is a header-only library that depends only on [`Eigen`](https://eigen.tuxfamily.org/index.php?title=Main_Page) (version >= 3.4.0), so the quick and dirty way of installing it is by simply copying the `include` directory to your project and make sure to have [`Eigen`](https://eigen.tuxfamily.org/index.php?title=Main_Page) available however you see fit. Alternatively, you can do things properly and use `CMake` (version >= 3.14).

Optionally, you can also install [`Sandals`](https://StoccoDavide.github.io/Sandals) ODEs/DAEs C++17 integrator library, and [`Matplot++`](https://alandefreitas.github.io/matplotplusplus/) (version >= 1.2.0) to plot the results of your simulations.

### CMake

If you are using CMake, you can add the library as a subdirectory in your project.

```cmake
add_subdirectory(path/to/Astro)
target_link_libraries(your_target PRIVATE Astro::Astro)
```

You can use `FetchContent` to download the library from GitHub.

```cmake
include(FetchContent)

# Optionally specify a custom path to fetch content to
set(FETCHCONTENT_BASE_DIR "path/to/your/dependencies")
fetchcontent_declare(
  Astro
  GIT_REPOSITORY https://github.com/StoccoDavide/Astro.git
  GIT_TAG        main
)
fetchcontent_makeavailable(Astro)
target_link_libraries(your_target PRIVATE Astro::Astro)
```

If you already have `Astro` somewhere on your system, you can use `find_pacakge` directly.

```cmake
# Optionally specify a custom path to find content from
list(APPEND CMAKE_PREFIX_PATH "path/to/your/dependencies")
find_package(
  Astro
  ${YOUR_DESIRED_ASTRO_VERSION}
  NO_MODULE
)

target_link_libraries(your_target PRIVATE Astro::Astro)
```

Since we are nice people, we also show you how to conditionally use `FetchContent` based if you already have the library or not.

```cmake
# Optionally specify a custom path to find content from
list(APPEND CMAKE_PREFIX_PATH "path/to/your/dependencies")
find_package(
  Astro
  ${YOUR_DESIRED_ASTRO_VERSION}
  NO_MODULE
)

if(NOT TARGET Astro::Astro)
  include(FetchContent)

  # Optionally specify a custom path to fetch content to
  set(FETCHCONTENT_BASE_DIR "path/to/your/dependencies")
  fetchcontent_declare(
    Astro
    GIT_REPOSITORY https://github.com/StoccoDavide/Astro.git
    GIT_TAG        main
  )

  fetchcontent_makeavailable(Astro)
endif()

target_link_libraries(your_target PRIVATE Astro::Astro)
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
