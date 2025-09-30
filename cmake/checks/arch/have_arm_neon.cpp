#include <arm_neon.h>

int main() {
  float32x4_t a = vdupq_n_f32(1.0f);
  return (int)vgetq_lane_f32(a, 0);
}
