#include <iostream>
#include <algorithm>
#include <limits>
#include <vector>
#include <set>

using namespace std;

template <typename T>
void print_mat(T **array, size_t rows, size_t cols)
{
    std::cout << __func__ << std::endl;
    for (size_t i = 0; i < rows; ++i)
    {
        std::cout << i << ": ";
        for (size_t j = 0; j < cols; ++j)
            std::cout << array[i][j] << '\t';
        std::cout << std::endl;
    }
}

vector<vector<int>> combination(const int N, const int M){
    vector<vector<int>> idx;
    
    if(N == M){
        vector<int> comb;
        for(int i = 0; i < M; ++i)
            comb.push_back(i);
        idx.push_back(comb);
        return idx;
    }
    
    vector<int> n, ind;

    for(int i = 0; i < M; ++i){
        n.emplace_back(i);
    }

    for(int i = 0; i < M - N; ++i)
        ind.emplace_back(0);
    for(int i = 0; i < N ; ++i)
        ind.emplace_back(1);

    // combination
    do{
        vector<int> comb;
        for(int i = 0; i < ind.size(); ++i){
            if(ind[i] == 1){
                comb.emplace_back(i);
            }
        }

        idx.emplace_back(comb);
    }while(next_permutation(ind.begin(), ind.end()));
    
    return idx;
}

template <size_t N, size_t M>
void createMat(const float (&src)[N][M], float ***arr, int& rows, int& cols, bool& bTrans, const int Mode){
    int mult = Mode == 0 ? 1 : -1;
    
    if(rows > cols){
        int temp = rows;
        rows = cols;
        cols = temp;
        bTrans = true;
    }
    
    *arr = new float *[rows];
    for(int i = 0; i < rows; ++i)
        (*arr)[i] = new float[cols];
    
    for(int i = 0; i < rows; ++i){
        for(int j = 0; j < cols; ++j){
            if(bTrans){
                (*arr)[i][j] = src[j][i] * mult;
            }
            else
                (*arr)[i][j] = src[i][j] * mult;
        }
    }
}

void init_solve(float ***arr, size_t rows, size_t cols){
    for(int i = 0; i < rows; ++i){
        float min = numeric_limits<float>::max();
        for(int j = 0; j < cols; ++j){
            if(min > (*arr)[i][j])
                min = (*arr)[i][j];
        }
        
        for(int j = 0; j < cols; ++j)
            (*arr)[i][j] -= min;
    }
    
    for(int i = 0; i < cols; ++i){
        float min = numeric_limits<float>::max();
        for(int j = 0; j < rows; ++j){
            if(min > (*arr)[j][i])
                min = (*arr)[j][i];
        }
        
        for(int j = 0; j < rows; ++j)
            (*arr)[j][i] -= min;
    }
}

bool dfs(float ***arr, size_t size, set<int>& colSet, int curRow, vector<pair<int, int>> &vecPair){
    for(int i = 0; i < size; ++i){
        if((*arr)[curRow][i] == 0 && colSet.find(i) == colSet.end()){
            vecPair.emplace_back(make_pair(curRow,i));
            colSet.insert(i);
            
            if(size-1 == curRow)
                return true;
            else if(dfs(arr, size, colSet, curRow+1, vecPair))
                return true;
            
            vecPair.pop_back();
            colSet.erase(i);
        }
    }
    
    return false;
}

bool isClear(float ***arr, size_t size, vector<pair<int, int>>& vecPair){
    //dfs
    set<int> colSet;
    
    vecPair.clear();
    if(dfs(arr, size, colSet, 0, vecPair))
        return true;
    else{
        //greedy allocation for fail case
        colSet.clear();
        vecPair.clear();
        for(int i = 0; i < size; ++i){
            for(int j = 0; j < size; ++j){
                if((*arr)[i][j] == 0 && colSet.find(j) == colSet.end()){
                    colSet.insert(j);
                    vecPair.emplace_back(make_pair(i,j));
                    break;
                }
            }
        }
    }
    
    return false;
}

// Now to the drawing part.
void refine_mat(float ***arr, size_t rows, size_t cols, vector<pair<int, int>>& pairs){
    set<int> uncoverdRows, coveredCols;
    vector<int> newlymarkedRows, newlymarkedCols;
    set<pair<int, int>> assignment;
    
    // Mark all rows having no assignments
    for(int i = 0, j = 0; i < rows; ++i){
        if(pairs.size() < j + 1
           || pairs.at(j).first != i)
            newlymarkedRows.push_back(i);
        else
            assignment.insert(pairs.at(j++));
    }
    
    while(true){
        if(newlymarkedRows.empty())
            break;
        
        // Mark all columns having zeros in newly marked row(s)
        newlymarkedCols.clear();
        for(auto& r:newlymarkedRows){
            uncoverdRows.insert(r);
            for(int i = 0; i < cols; ++i){
                if((*arr)[r][i] == 0
                   && coveredCols.find(i) == coveredCols.end())
                    newlymarkedCols.push_back(i);
            }
        }
        
        if(newlymarkedCols.empty())
            break;
        
        // Mark all rows having assignments in newly marked columns
        newlymarkedRows.clear();
        for(auto& c:newlymarkedCols){
            coveredCols.insert(c);
            for(int i = 0; i < rows; ++i){
                if((*arr)[i][c] == 0
                   && assignment.find(make_pair(i,c)) != assignment.end()
                   && uncoverdRows.find(i) == uncoverdRows.end())
                    newlymarkedRows.push_back(i);
            }
        }
    }
    
    // Now draw lines through all marked columns and unmarked rows.
    // From the elements that are left, find the lowest value.
    float minval = numeric_limits<float>::max();
    for(int i = 0; i < rows; ++i){
        for(int j = 0; j < cols; ++j){
            if(uncoverdRows.find(i) != uncoverdRows.end()
               && coveredCols.find(j) == coveredCols.end()
               && minval > (*arr)[i][j])
                minval = (*arr)[i][j];
        }
    }
    
    // Subtract this from every unmarked element and add it to every element covered by two lines.
    for(int i = 0; i < rows; ++i){
        for(int j = 0; j < cols; ++j){
            if(uncoverdRows.find(i) != uncoverdRows.end()
               && coveredCols.find(j) == coveredCols.end())
                (*arr)[i][j] -= minval;
            
            if(uncoverdRows.find(i) == uncoverdRows.end()
               && coveredCols.find(j) != coveredCols.end())
               (*arr)[i][j] += minval;
        }
    }
    
}

float partial_solve(float ***array, size_t rows, size_t cols, vector<int>* src_idx, int Mode, vector<int>& assignment_idx){
    float **arr = nullptr;
    vector<pair<int, int>> pairs;
    
    arr = new float *[rows];
    for(int i = 0; i < rows; ++i){
        arr[i] = new float[src_idx->size()];
        for(int j = 0; j < src_idx->size(); ++j)
            arr[i][j] = (*array)[i][src_idx->at(j)];
    }
    
    print_mat(arr, src_idx->size(), src_idx->size());
    init_solve(&arr, src_idx->size(), src_idx->size());
    print_mat(arr, src_idx->size(), src_idx->size());
    
    while(true){
        if(isClear(&arr, rows, pairs))
            break;
        
        refine_mat(&arr, rows, cols, pairs);
        print_mat(arr, src_idx->size(), src_idx->size());
    }
    
    for (int i = 0; i < rows; ++i)
        delete[] arr[i];
    delete[] arr;
    
    //finalize
    float val = 0;
    int tempidx = 0;
    for(auto& p:pairs)
        assignment_idx.push_back(src_idx->at(p.second));
    for(auto& idx:assignment_idx)
        val += (*array)[tempidx++][idx];
    
    return val;
}

template <size_t N, size_t M>
float Solve(const float (&Cost)[N][M], const int MODE, vector<int>& assignment_index){
    float val = numeric_limits<float>::max();
    assignment_index.clear();
    assignment_index.reserve(N);
    
    float **arr = nullptr;
    int rows = N;
    int cols = M;
    bool bTrans = false;
    
    createMat(Cost, &arr, rows, cols, bTrans, MODE);
    print_mat(arr, rows, cols);
    
    vector<vector<int>> comb = combination(rows, cols);
    for(auto& c:comb){
        vector<int> partial_assign;
        float partial_val = partial_solve(&arr, rows, cols, &c, MODE, partial_assign);
        
        if(partial_val < val){
            val = partial_val;
            assignment_index = partial_assign;
        }
    }
    
    if(bTrans){
        int *transposed = new int [N];
        int cnt = 0;
        fill(transposed, transposed + N, -1);
        for(auto& idx:assignment_index)
            transposed[idx] = cnt++;
        
        assignment_index.clear();
        for(int i = 0; i < N; ++i)
            assignment_index.push_back(transposed[i]);
        
        delete[] transposed;
    }
    
    return MODE ? val * -1 : val;
}

int main(int argc, const char * argv[]) {
    int N = 3;
    int M = 4;
    vector<int> assignment_idx;
    
//    float Cost[3][4] =
//    {
//        {3,7,5,11},
//        {5,4,6,3},
//        {6,10,1,1}
//    };
    
//    float Cost[3][3] =
//    {
//        {3, 8, 9},
//        {4, 12, 7},
//        {4, 8, 5},
//    };
    
    float Cost[4][3] = {
        {3, 5, 6},
        {7, 4, 10},
        {5, 6, 1},
        {11, 3, 1}
    };
    
    float optCost = Solve(Cost, 0, assignment_idx);
    cout << "Total Cost : " << optCost << ", assignment_index = { ";
    for(auto& idx:assignment_idx){
        cout << idx << " ";
    }
    cout << "}" << endl;
    
    return 0;
}
