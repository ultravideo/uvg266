#include <altivec.h>

int main() {
  vector unsigned int a = (vector unsigned int){1, 2, 3, 4};
  vector unsigned int b = (vector unsigned int){5, 6, 7, 8};
  vector unsigned int c = vec_add(a, b);
  return vec_extract(c, 0);
}
