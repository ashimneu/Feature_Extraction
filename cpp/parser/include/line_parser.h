#ifndef ECS_CLIENT_LINE_PARSER_H
#define ECS_CLIENT_LINE_PARSER_H

size_t
parse_block(const char * text_line, size_t start_pos, char * sub_string, size_t max_block_size, char delimiter)
{
    // reset sub_string
    memset(sub_string, 0, max_block_size);
    size_t j {0};
    size_t pos;
    for(pos = start_pos; pos < start_pos + max_block_size; ++pos)
    {
        if(*(text_line + pos) == delimiter)
        {
            return ++pos;
        }
        else
        {
            *(sub_string+j) += *(text_line + pos);
            ++j;
        }
    }
    return pos;
}

#endif // ECS_CLIENT_LINE_PARSER_H
