#include <nmmintrin.h>

int main() {
  __m128i a = _mm_set1_epi8(0);
  __m128i b = _mm_set1_epi8(0);
  return _mm_cmpestri(a, 0, b, 0, _SIDD_CMP_EQUAL_ANY);
}
