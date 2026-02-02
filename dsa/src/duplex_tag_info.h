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

#endif
