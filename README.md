# Matrix

C++17 header which simplifies the use of matrices.

## Usage examples

#### Including the header file

```cpp
#include "matrix.h"
```

### Basic

#### Creating

+ a matrix filled with your elements

```cpp
nur::Matrix<int> m{
    {
        {
            4, 4,
        },
        {
            6, 5,
        },
    }
};
```

+ a matrix of size `3x4`

```cpp
nur::Matrix<int> m(3, 4);
```

#### Element access using `operator []`

```cpp
nur::Matrix<double> m(4, 5);
for (int i = 2; i < 4; ++i) {
    for (it j = 1; j < 3; ++j) {
        m[i][j] = 3.14;
    }
}
```
//TODO
