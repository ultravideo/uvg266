#include <catch2/catch_session.hpp>

auto main(int argc, char** argv) -> int
{
  Catch::Session session;
  if (session.applyCommandLine(argc, argv) != 0) {
    return EXIT_FAILURE;
  }

  return session.run();
}
