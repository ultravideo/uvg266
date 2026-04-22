#include <emmintrin.h>

int main() {
  __m128i a = _mm_setzero_si128();
  __m128i b = _mm_set1_epi32(1);
  __m128i c = _mm_add_epi32(a, b);
  return _mm_cvtsi128_si32(c);
}
