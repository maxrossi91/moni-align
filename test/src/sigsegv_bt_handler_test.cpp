/* sigsegv_bt_handler_test - Test bt_handler catches SIGSEGV errors
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
   \file sigsegv_bt_handler_test.cpp
   \brief sigsegv_bt_handler_test.cpp Test bt_handler catches SIGSEGV errors.
   \author Massimiliano Rossi
   \date 17/11/2020
*/

#include <iostream>

#define VERBOSE

#include <common.hpp>
#include <stacktrace.hpp>
#include <signal.h>

int func_a(int a, char b) {

  char *p = (char *)0xdeadbeef;

  a = a + b;
  *p = 10;  /* CRASH here!! */

  return 2*a;
}


int sigsegv_trigger() {

  int res, a = 5;

  res = 5 + func_a(a, 't');

  return res;
}

int main(int argc, char *const argv[])
{
  enable_stacktrace();
  sigsegv_trigger();
  return 0;
}