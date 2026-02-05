#ifndef DUPLEX_TAG_INFO_H_
#define DUPLEX_TAG_INFO_H_

#include <sstream>
#include <string>
#include <vector>

typedef struct duplex_tag_info_t {
    int32_t beg;
    int32_t end;
    std::string fwd_bc;
    std::string rev_bc;
} duplex_tag_info_t;

static bool duplex_tag_info_is_pos_in_template(const duplex_tag_info_t *info, const int32_t pos) {
    // ASSUMPTION: offset correction (based on the convention used for the position)
    //  has been applied to the input position.
    return (pos >= info->beg) && (pos <= info->end);
}

static duplex_tag_info_t duplex_tag_info_parse(std::string idf1) {
    std::istringstream iss(idf1);
    std::vector<std::string> tokens;
    std::string token;

    // TODO: optimise!
    while (std::getline(iss, token, ',')) {
        tokens.push_back(token);
    }
    duplex_tag_info_t idf = {
        .beg = std::stoi(tokens[1]),
        .end = std::stoi(tokens[2]),
        .fwd_bc = tokens[3],
        .rev_bc = tokens[4]
    };
    return idf;
}

#endif
