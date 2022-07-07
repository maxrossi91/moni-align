/* backward_test - Test of the backward library
    Copyright (C) 2020 Massimiliano Rossi
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*!
   \file backward_test.cpp
   \brief backward_test.cpp Test of the backward library.
   \author Massimiliano Rossi
   \date 17/11/2020
*/

#include <backward.hpp>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace backward;

class TracedException : public std::runtime_error {
public:
  TracedException() : std::runtime_error(_get_trace()) {}

private:
  std::string _get_trace() {
    std::ostringstream ss;

    StackTrace stackTrace;
    TraceResolver resolver;
    stackTrace.load_here();
    resolver.load_stacktrace(stackTrace);

    for (std::size_t i = 0; i < stackTrace.size(); ++i) {
      const ResolvedTrace trace = resolver.resolve(stackTrace[i]);

      ss << "#" << i << " at " << trace.object_function << "\n";
    }

    return ss.str();
  }
};

void f(int i) {
  if (i >= 42) {
    throw TracedException();
  } else {
    std::cout << "i=" << i << "\n";
    f(i + 1);
  }
}

int main() {
  try {
    f(0);
  } catch (const TracedException &ex) {
    std::cout << ex.what();
  }
}