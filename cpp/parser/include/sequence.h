#ifndef ECS_CLIENT_SEQUENCE_H
#define ECS_CLIENT_SEQUENCE_H

std::vector<int> sequence(int start, int end, int increment){
    std::vector<int> sequence;

    if (increment > 0 ){
        for (int i = start; i <= end; i +=increment){
            sequence.push_back(i);
        }
    }
    else{
        for (int i = start; i >= end; i += increment){
            sequence.push_back(i);
        }
    }
    return sequence;
}
#endif // ECS_CLIENT_SEQUENCE_H