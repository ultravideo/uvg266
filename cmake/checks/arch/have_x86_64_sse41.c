#include <smmintrin.h>

int main() {
  __m128i a = _mm_setzero_si128();
  a = _mm_insert_epi32(a, 1, 0);
  return _mm_cvtsi128_si32(a);
}
