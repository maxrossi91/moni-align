/* stacktrace - Handle stacktrace print on accidental crashes
    Copyright (C) 2022 Massimiliano Rossi

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
   \file stacktrace.hpp
   \brief stacktrace.hpp Handle stacktrace print on accidental crashes.
   \author Massimiliano Rossi
   \date 05/05/2022
*/

#ifndef _STACKTRACE_HH
#define _STACKTRACE_HH

// #define BACKWARD_HAS_DW 1
// #define BACKWARD_HAS_BFD 1

#include <iostream>
#include <execinfo.h>
#include <cxxabi.h>
#include <signal.h>
#include <string>
#include <regex>
#include <ucontext.h>
#include <backward.hpp>

// //*********************** Stack trace options *******************************
// Inspired from https://panthema.net/2008/0901-stacktrace-demangled/

// Other sources: https://stackoverflow.com/questions/3355683/c-stack-trace-from-unhandled-exception
// Other sources: https://stackoverflow.com/questions/4636456/how-to-get-a-stack-trace-for-c-using-gcc-with-line-number-information
// Other sources: https://stackoverflow.com/questions/691719/c-display-stack-trace-on-exception

typedef struct trace_line_t{
    std::string file_name = "";
    std::string function_name = "";
    std::string offset = "";
    std::string address = "";

    friend std::ostream& operator<<(std::ostream& os, const struct trace_line_t& line);
} trace_line_t;

std::ostream& operator<<(std::ostream& os, const struct trace_line_t& line)
{
    os << "\t" << line.address << "\t" << line.file_name << " : " << line.function_name << "+" << line.offset;
    return os;
}


static inline trace_line_t parse_symbol(char* symbol)
{

    trace_line_t line;
	// find parentheses and +address offset surrounding the mangled name:
	// ./module(function+0x15c) [0x8048a6d]

    char* start = symbol;
    int state = 0;
	for (char *p = symbol; *p; ++p)
	{
	    if (*p == '(' && state == 0) 
            *p = '\0', line.file_name = std::string(start), start = p+1, state++;
        else if (*p == '+' && state == 1)
            *p = '\0', line.function_name = std::string(start), start = p + 1, state++;
	    else if (*p == ')' && state == 2) 
            *p = '\0', line.offset = std::string(start), start = p + 1, state++;
	}
    line.address = std::string(start + 1); // +1 to skipp the space

    if(state == 3)
    {
        // allocate string which will be filled with the demangled function name
        size_t funcnamesize = 256;
        char *funcname = (char *)malloc(funcnamesize);

        // mangled name is now in [begin_name, begin_offset) and caller
        // offset in [begin_offset, end_offset). now apply
        // __cxa_demangle():

        int status;
        char *ret = abi::__cxa_demangle(line.function_name.c_str(), funcname, &funcnamesize, &status);

        if (status == 0) line.function_name = std::string(ret);
        else line.function_name += "()";

        free(funcname);
    }

    return line;
}

static inline void print_stacktrace(std::ostream& stream = std::cerr, unsigned int max_frames = 63)
{
    stream << "Stack trace:" << std::endl;

    // storage array for stack trace address data
    void* addrlist[max_frames+1];

    // retrieve current stack addresses
    int addrlen = backtrace(addrlist, sizeof(addrlist) / sizeof(void*));

    if (addrlen == 0) {
        stream <<  "\t<empty, possibly corrupt>\n" << std::endl;
    	return;
    }

    // resolve addresses into strings containing "filename(function+address)",
    // this array must be free()-ed
    char** symbollist = backtrace_symbols(addrlist, addrlen);



    // iterate over the returned symbol lines. skip the first, it is the
    // address of this function.
    for (int i = 1; i < addrlen; i++)
        stream << parse_symbol(symbollist[i]) << std::endl;
	
    free(symbollist);
}

void stacktrace_sighandler(int sig)
{
  print_stacktrace();
  exit(sig);
}

// Other option using backward-cpp library
void stacktrace_sigaction(int signalNumber, siginfo_t *signalInfo, void *signalContext) {
    // This holds the context that the signal came from, including registers and stuff
    ucontext_t* context = (ucontext_t*) signalContext;

    // Fetch out the registers
    void* ip = (void*)context->uc_mcontext.gregs[REG_RIP];
    void** bp = (void**)context->uc_mcontext.gregs[REG_RBP];
    
    static backward::StackTrace stack_trace;

    stack_trace.load_from(ip, 32);

    static backward::Printer p;
    p.object = true;
    p.color_mode = backward::ColorMode::automatic;
    p.address = true;
    p.print(stack_trace, std::cerr);

    exit(signalNumber);
}


static inline void enable_stacktrace(){
    struct sigaction sa;

    memset(&sa, 0, sizeof(struct sigaction));
    sigemptyset(&sa.sa_mask);
    // sa.sa_handler = stacktrace_sighandler;
    sa.sa_sigaction = stacktrace_sigaction;
    sa.sa_flags = SA_SIGINFO;

    sigaction(SIGABRT, &sa, nullptr);
    sigaction(SIGSEGV, &sa, nullptr);
    sigaction(SIGBUS, &sa, nullptr);
    sigaction(SIGILL, &sa, nullptr);
    sigaction(SIGFPE, &sa, nullptr);

}

#endif /* end of include guard: _STACKTRACE_HH */
