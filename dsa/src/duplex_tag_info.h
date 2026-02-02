#ifndef DUPLEX_TAG_INFO_H_
#define DUPLEX_TAG_INFO_H_

#include <sstream>
#include <string>

typedef struct duplex_tag_info {
    int beg;
    int end;
    std::string fwd_bc;
    std::string rev_bc;
} duplex_tag_info;

static duplex_tag_info duplex_tag_info_parse(std::string idf1) {
    std::istringstream iss(idf1);
    std::vector<std::string> tokens;
    std::string token;

    // TODO: optimise!
    while (std::getline(iss, token, ',')) {
        tokens.push_back(token);
    }
    duplex_tag_info idf = {
        .beg = std::stoi(tokens[1]),
        .end = std::stoi(tokens[2]),
        .fwd_bc = tokens[3],
        .rev_bc = tokens[4]
    };
    return idf;
}

typedef struct duplex_info_t {
    // std::string id;
    uint64_t index;  // use as key
    uint8_t type = 0;  // filter out on zero
} duplex_info_t;

/*
std::map<std::string, duplex_info_t> duplex_id_infos;
std::unsorted_map<uint64_t, duplex_info_t> decoder;
foreach k, v:
    if type != 0:
        decoder[v.index] = v
*/

#endif
