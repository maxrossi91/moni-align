#include <vector>
#include <string>
#include <fstream>

//logger initialized per execution, each thread should have access, where to put???
struct logger_t {
    std::string filename = ; 
    std::vector<std::string> reads;

    void write() {
        std::fstream file;
        file.open(filename, std::ios::out);
        file << "[ ";
        for(int i = 0; i < reads.size(); i++)
        {
            file << " {" << std::endl;
            file << reads[i];
            if(i == reads.size()-1) {
                file <<"} ";
            }
            else {
                file <<"},";
            }
        }
        file << "]";

    }
};

//This is created and updated for every read and then added to logger
//what information does max need for the debugger?
//what is kseq_t