#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include "beam_cky_parser.h"
#include "beam_cky_parser.cc"
#include "../Utils/reader.h"
#include "../Utils/common.h"
#include "../Utils/codon.h"

// #ifndef CODON_TABLE
// #define CODON_TABLE "./codon_usage_freq_table_human.csv"
// #endif

#ifndef CODING_WHEEL
#define CODING_WHEEL "./coding_wheel.txt"
#endif

#ifndef GAP_PENALTY
#define GAP_PENALTY 5
#endif

#ifndef STEM_THRESHOLD
#define STEM_THRESHOLD 30
#endif
// also should fix in forgi/forgi_run.py

using namespace LinearDesign;

template <typename ScoreType, typename IndexType>
bool output_result(const DecoderResult<ScoreType, IndexType>& result, 
        const double duration, const double lambda, const bool is_verbose, 
        const Codon& codon, string& CODON_TABLE, string& utr3_seq, const string& penalties) {

    stringstream ss;
    string result_cds = result.sequence;
    string result_utr;
    string result_seq = result.sequence;
    
    if (!utr3_seq.empty()) {
        result_cds = result.sequence.substr(0, result.sequence.length() - utr3_seq.length());
        std::transform(utr3_seq.begin(), utr3_seq.end(), utr3_seq.begin(), ::tolower);
        result_seq = result_cds + utr3_seq;
    }

    if (is_verbose)
        ss << "Using lambda = " << (lambda / 100.) << "; Using codon frequency table = " << CODON_TABLE << endl;
    ss << "mRNA sequence:  " << result_seq << endl;
    ss << "mRNA structure: " << result.structure << endl;
    ss << "mRNA folding free energy: " << std::setprecision(2) << fixed << result.score << " kcal/mol" << endl;
    ss << "mRNA CAI: " << std::setprecision(3) << fixed << codon.calc_cai(result_cds) << endl;
    if (is_verbose)
        ss << "Runtime: " << duration << " seconds" << endl;
    if(!penalties.empty()){
        cout << "Penalties: " << penalties << endl;
    }
    cout << ss.str() << endl;

    return true;
}

void show_usage() {
    cerr << "echo SEQUENCE | ./lineardesign -l [LAMBDA] -v -c [codon table] -u3 [3' UTR] -p [penalty string] -g [gap to search]" << endl;
    cerr << "OR" << endl;
    cerr << "cat SEQ_FILE_OR_FASTA_FILE | ./lineardesign -l [LAMBDA] -v -c [codon table] -u3 [3' UTR] -p [penalty string] -g [gap to search]" << endl;
}

//==========================================================================
BeamCKYParser<ScoreType, IndexType>& parse_penalty_options(BeamCKYParser<ScoreType, IndexType>& parser, string penalty_str) {
    if (penalty_str.empty())
        return parser;
    parser.penalty_mode = PENALTY_ON;

    std::istringstream stream(penalty_str);
    std::string pair_str;

    while (std::getline(stream, pair_str, ',')){
        size_t tilde_pos = pair_str.find('~');
        if(tilde_pos != std::string::npos){
            IndexType penalty_start = std::stoi(pair_str.substr(0, tilde_pos));
            IndexType penalty_end = std::stoi(pair_str.substr(tilde_pos + 1));
            if(penalty_start > penalty_end){
                IndexType temp = penalty_end;
                penalty_end = penalty_start;
                penalty_start = temp;
            }
            parser.penalty_vec.emplace_back(penalty_start, penalty_end);
        }
        else {
            std::cerr << "Invalid format in penalty_str: " << pair_str << endl;
            exit(EXIT_FAILURE);
        }
    }

    return parser;
}

//==========================================================================

int main(int argc, char** argv) {
    // default args
    double lambda = 0.0f;
    bool is_verbose = false;
    string CODON_TABLE = "./codon_usage_freq_table_human.csv";
    string utr3_seq = "";
    int starnum = 0;

    // parse args
    if (argc != 6) {
        show_usage();
        return 1;
    }else{
        lambda = atof(argv[1]);
        is_verbose = atoi(argv[2]) == 1;
        if (string(argv[3]) != ""){
            CODON_TABLE = argv[3];
        }
        utr3_seq = string(argv[4]);
        utr3_seq = utr3_seq.substr(0, utr3_seq.length() - (utr3_seq.length() % 3));
    } 
    lambda *= 100.;
    
    // load codon table and coding wheel
    Codon codon(CODON_TABLE);
    std::unordered_map<string, Lattice<IndexType>> aa_graphs_with_ln_weights;
    std::unordered_map<std::string, std::unordered_map<std::tuple<NodeType, NodeType>, std::tuple<double, NucType, NucType>, std::hash<std::tuple<NodeType, NodeType>>>> best_path_in_one_codon_unit;
    std::unordered_map<std::string, std::string> aa_best_path_in_a_whole_codon;
    prepare_codon_unit_lattice<IndexType>(CODING_WHEEL, codon, aa_graphs_with_ln_weights, best_path_in_one_codon_unit, aa_best_path_in_a_whole_codon, lambda);

    // main loop
    string aa_seq, aa_tri_seq;
    vector<string> aa_seq_list, aa_name_list;
    // load input
    for (string seq; getline(cin, seq);){
        if (seq.empty()) continue;
        if (seq[0] == '>'){
            aa_name_list.push_back(seq); // sequence name
            if (!aa_seq.empty())
                aa_seq_list.push_back(aa_seq);
            aa_seq.clear();
            continue;
        }else{
            rtrim(seq);
            aa_seq += seq;
        }
    }
    if (!aa_seq.empty())
        aa_seq_list.push_back(aa_seq);

    // start design
    for(int i = 0; i < aa_seq_list.size(); i++){
        auto& aa_seq = aa_seq_list[i];
        if (!utr3_seq.empty())
            aa_seq = aa_seq + codon.cvt_rna_seq_to_aa_seq(utr3_seq);
        // convert to uppercase
        transform(aa_seq.begin(), aa_seq.end(), aa_seq.begin(), ::toupper);
        aa_tri_seq.clear();
        if (is_verbose)
            cout << "Input protein: " << aa_seq << endl;
        if (!ReaderTraits<Fasta>::cvt_to_seq(aa_seq, aa_tri_seq)) 
            continue;

        // init parser
        BeamCKYParser<ScoreType, IndexType> parser(lambda, is_verbose);
        

        auto protein = util::split(aa_tri_seq, ' ');
        // parse
        auto system_start = chrono::system_clock::now();
        
        //==========================================================================
        DFA<IndexType> dfa;
        if (utr3_seq.empty()) {
            dfa = get_dfa<IndexType>(aa_graphs_with_ln_weights, util::split(aa_tri_seq, ' '));
        } else {
            string command = "python3 src/dfa/dfa_gen.py " + aa_seq + " -l " + std::to_string(lambda) + " -u3 " + utr3_seq;
            int dfa_result = system(command.c_str());
            dfa = get_dfa<IndexType>("./src/dfa/dfa.txt", utr3_seq);
        }
        parser = parse_penalty_options(parser, string(argv[5]));
        //==========================================================================

        auto result = parser.parse(dfa, codon, aa_seq, protein, aa_best_path_in_a_whole_codon, best_path_in_one_codon_unit, aa_graphs_with_ln_weights);
        auto system_diff = chrono::system_clock::now() - system_start;
        auto system_duration = chrono::duration<double>(system_diff).count();
        output_result(result, system_duration, lambda, is_verbose, codon, CODON_TABLE, utr3_seq, string(argv[5]));


#ifdef FINAL_CHECK
        if (codon.cvt_rna_seq_to_aa_seq(result.sequence) != aa_seq) {
            std::cerr << "Final Check Failed:" << std::endl;
            std::cerr << codon.cvt_rna_seq_to_aa_seq(result.sequence) << std::endl;
            std::cerr << aa_seq << std::endl;
            assert(false);
        }
#endif
    }
    return 0;
}