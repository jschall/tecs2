#include "DFParser.h"
#include <iostream>
#include "json.hpp"
#include <fcntl.h>
#include <sys/mman.h>
#include <iomanip>
#include <sstream>
#include <boost/algorithm/string/join.hpp> // Include for boost::split
#include <math.h>

using namespace std;
using namespace nlohmann;

int main(int argc, char** argv) {
    if (argc < 2) {
        cout << "too few arguments" << endl;
        return 1;
    }

    std::ios::sync_with_stdio(true);

    int fp = open(argv[1], O_RDONLY, 0);
    size_t logdata_len = lseek(fp, 0, SEEK_END);
    lseek(fp, 0, SEEK_SET);
    uint8_t* logdata = (uint8_t*)mmap(0, logdata_len, PROT_READ, MAP_SHARED, fp, 0);
    madvise(logdata, logdata_len, POSIX_MADV_SEQUENTIAL);

    DFParser parser(logdata, logdata_len);

    uint64_t count=0;
    uint8_t type;

    float rho=0;
    float TAS=0;
    float prop_omega=0;
    float STEdot=0;
    float bank_rad=0;

    bool have_TECS = false;
    bool have_CTUN = false;
    bool have_ATT = false;
    bool have_ESC = false;

    cout << "rho,TAS,prop_omega,STEdot,bank_rad" << endl;
    DFParser::message_t msg;
    while(parser.next_message(msg)) {
        auto& fields = parser.get_fields(msg);

        if (parser.get_message_name(msg) == "ESC") {
            uint32_t instance;
            parser.get_scalar_field(msg,"Instance", instance);
            if (instance == 4) {
                float rpm;
                parser.get_scalar_field(msg,"RPM", rpm);
                prop_omega = rpm*2*M_PI/60;
                have_ESC = true;
            }
        }

        if (parser.get_message_name(msg) == "TECS") {
            float sp;
            parser.get_scalar_field(msg,"sp", sp);
            TAS = sp;
            have_TECS = true;
        }

        if (parser.get_message_name(msg) == "CTUN") {
            float E2T;
            parser.get_scalar_field(msg,"E2T", E2T);
            rho = 1.225/sqrt(E2T);
            have_CTUN = true;
        }

        if (parser.get_message_name(msg) == "ATT") {
            float roll;
            parser.get_scalar_field(msg,"roll", roll);
            bank_rad = M_PI/180.0 * roll;
            have_ATT = true;
        }

        if (parser.get_message_name(msg) == "TEC3") {
            float KED;
            float PED;
            parser.get_scalar_field(msg,"KED", KED);
            parser.get_scalar_field(msg,"PED", PED);
            STEdot = KED+PED;

            if (have_TECS && have_CTUN && have_ATT && have_ESC) {
                cout << rho << "," << TAS << "," << prop_omega << "," << STEdot << "," << bank_rad << endl;
            }
        }

//         cout << parser.get_message_name(msg) << " {";
//
//         bool first_iter = true;
//         for (auto field : fields) {
//             if (!first_iter) {
//                 cout << ", ";
//             } else {
//                 first_iter = false;
//             }
//
//             cout << field.name << " : " << parser.get_value_string(msg,field);
//         }
//
//         cout << "}" << endl;
    }

    return 0;
}
