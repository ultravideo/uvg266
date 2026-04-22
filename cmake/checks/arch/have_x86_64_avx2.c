#include <immintrin.h>

int main() {
  __m256i a = _mm256_setzero_si256();
  __m256i b = _mm256_set1_epi32(1);
  __m256i c = _mm256_add_epi32(a, b);
  return _mm256_extract_epi32(c, 0);
}
