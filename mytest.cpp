// CMSC 341 - Spring 2025 - Project 4
#include "dnadb.h"
#include <math.h>
#include <algorithm>
#include <random>
#include <vector>
using namespace std;
enum RANDOM {UNIFORMINT, UNIFORMREAL, NORMAL, SHUFFLE};
class Random {
public:
    Random(){}
    Random(int min, int max, RANDOM type=UNIFORMINT, int mean=50, int stdev=20) : m_min(min), m_max(max), m_type(type)
    {
        if (type == NORMAL){
            //the case of NORMAL to generate integer numbers with normal distribution
            m_generator = std::mt19937(m_device());
            //the data set will have the mean of 50 (default) and standard deviation of 20 (default)
            //the mean and standard deviation can change by passing new values to constructor 
            m_normdist = std::normal_distribution<>(mean,stdev);
        }
        else if (type == UNIFORMINT) {
            //the case of UNIFORMINT to generate integer numbers
            // Using a fixed seed value generates always the same sequence
            // of pseudorandom numbers, e.g. reproducing scientific experiments
            // here it helps us with testing since the same sequence repeats
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_unidist = std::uniform_int_distribution<>(min,max);
        }
        else if (type == UNIFORMREAL) { //the case of UNIFORMREAL to generate real numbers
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_uniReal = std::uniform_real_distribution<double>((double)min,(double)max);
        }
        else { //the case of SHUFFLE to generate every number only once
            m_generator = std::mt19937(m_device());
        }
    }
    void setSeed(int seedNum){
        // we have set a default value for seed in constructor
        // we can change the seed by calling this function after constructor call
        // this gives us more randomness
        m_generator = std::mt19937(seedNum);
    }
    void init(int min, int max){
        m_min = min;
        m_max = max;
        m_type = UNIFORMINT;
        m_generator = std::mt19937(10);// 10 is the fixed seed value
        m_unidist = std::uniform_int_distribution<>(min,max);
    }
    void getShuffle(vector<int> & array){
        // this function provides a list of all values between min and max
        // in a random order, this function guarantees the uniqueness
        // of every value in the list
        // the user program creates the vector param and passes here
        // here we populate the vector using m_min and m_max
        for (int i = m_min; i<=m_max; i++){
            array.push_back(i);
        }
        shuffle(array.begin(),array.end(),m_generator);
    }

    void getShuffle(int array[]){
        // this function provides a list of all values between min and max
        // in a random order, this function guarantees the uniqueness
        // of every value in the list
        // the param array must be of the size (m_max-m_min+1)
        // the user program creates the array and pass it here
        vector<int> temp;
        for (int i = m_min; i<=m_max; i++){
            temp.push_back(i);
        }
        std::shuffle(temp.begin(), temp.end(), m_generator);
        vector<int>::iterator it;
        int i = 0;
        for (it=temp.begin(); it != temp.end(); it++){
            array[i] = *it;
            i++;
        }
    }

    int getRandNum(){
        // this function returns integer numbers
        // the object must have been initialized to generate integers
        int result = 0;
        if(m_type == NORMAL){
            //returns a random number in a set with normal distribution
            //we limit random numbers by the min and max values
            result = m_min - 1;
            while(result < m_min || result > m_max)
                result = m_normdist(m_generator);
        }
        else if (m_type == UNIFORMINT){
            //this will generate a random number between min and max values
            result = m_unidist(m_generator);
        }
        return result;
    }

    double getRealRandNum(){
        // this function returns real numbers
        // the object must have been initialized to generate real numbers
        double result = m_uniReal(m_generator);
        // a trick to return numbers only with two deciaml points
        // for example if result is 15.0378, function returns 15.03
        // to round up we can use ceil function instead of floor
        result = std::floor(result*100.0)/100.0;
        return result;
    }

    string getRandString(int size){
        // the parameter size specifies the length of string we ask for
        // to use ASCII char the number range in constructor must be set to 97 - 122
        // and the Random type must be UNIFORMINT (it is default in constructor)
        string output = "";
        for (int i=0;i<size;i++){
            output = output + (char)getRandNum();
        }
        return output;
    }
    
    int getMin(){return m_min;}
    int getMax(){return m_max;}
    private:
    int m_min;
    int m_max;
    RANDOM m_type;
    std::random_device m_device;
    std::mt19937 m_generator;
    std::normal_distribution<> m_normdist;//normal distribution
    std::uniform_int_distribution<> m_unidist;//integer uniform distribution
    std::uniform_real_distribution<double> m_uniReal;//real uniform distribution
};

class Tester{
    public:
    bool linearInsertionTest();
    bool getDNAErrorCase();
    bool getDNAstandardcase();
    bool getDNAcollisioncase();
    bool removeTestStandardCase();
    bool rehashTest1();
    bool rehashTest2();
    bool rehashTest3();
    bool rehashTest4();
    bool quadraticInsertionTest();
    bool doublehashInsertionTest();
    bool removeErrorCase();
    bool updateLocNormalCase();
    bool updateLocErrorCase();
};

unsigned int hashCode(const string str);
string sequencer(int size, int seedNum);

int main(){
    Tester tester;
    cout << "Test 1: Linear Insertion Test" << endl;
    if(tester.linearInsertionTest()){
        cout << "test passed" << endl;
    }else{
        cout << "test failed" << endl;
    }
    cout << "Test 2: GetDNA error case Test" << endl;
    if(tester.getDNAErrorCase()){
        cout << "test passed" << endl;
    }else{
        cout << "test failed" << endl;
    }
    cout << "Test 3: GetDNA Test: no collisions" << endl;
    if(tester.getDNAstandardcase()){
        cout << "test passed" << endl;
    }else{
        cout << "test failed" << endl;
    }
    cout << "Test 4: GetDNA Test: multiple collisions" << endl;
    if(tester.getDNAcollisioncase()){
        cout << "test passed" << endl;
    }else{
        cout << "test failed" << endl;
    }    
    cout << "Test 5: Remove Test: No Collisions" << endl;
    if(tester.removeTestStandardCase()){
        cout << "test passed" << endl;
    }else{
        cout << "test failed" << endl;
    }
    cout << "Test 6: Rehash Test: First Rehash Iteration" << endl;
    if(tester.rehashTest1()){
        cout << "test passed" << endl;
    }else{
        cout << "test failed" << endl;
    }
    cout << "Test 7: Rehash Test: Rehashing Completion" << endl;
    if(tester.rehashTest2()){
        cout << "test passed" << endl;
    }else{
        cout << "test failed" << endl;
    }
    cout << "Test 8: Rehash Test: Rehashing After Deletion" << endl;
    if(tester.rehashTest3()){
        cout << "test passed" << endl;
    }else{
        cout << "test failed" << endl;
    }
    cout << "Test 9: Rehash Test: Rehashing Completely After Deletion" << endl;
    if(tester.rehashTest4()){
        cout << "test passed" << endl;
    }else{
        cout << "test failed" << endl;
    }
    cout << "Test 10: Quadratic Insertion Test" << endl;
    if(tester.quadraticInsertionTest()){
        cout << "test passed" << endl;
    }else{
        cout << "test failed" << endl;
    }
    cout << "Test 11: Double Hash Insertion Test" << endl;
    if(tester.doublehashInsertionTest()){
        cout << "test passed" << endl;
    }else{
        cout << "test failed" << endl;
    }
    cout << "Test 12: Remove Error Test" << endl;
    if(tester.removeErrorCase()){
        cout << "test passed" << endl;
    }else{
        cout << "test failed" << endl;
    }
    cout << "Test 13: Update Location Normal Case" << endl;
    if(tester.updateLocNormalCase()){
        cout << "test passed" << endl;
    }else{
        cout << "test failed" << endl;
    }
    cout << "Test 14: Update Location Error Case" << endl;
    if(tester.updateLocErrorCase()){
        cout << "test passed" << endl;
    }else{
        cout << "test failed" << endl;
    }
    return 0;
}
//tests insertion function//
bool Tester::linearInsertionTest(){
    Random RndLocation(MINLOCID,MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode, LINEAR);
    for (int i=0;i<10;i++){
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum(), true);
        dnadb.insert(dataObj);
    }
    //adding data points with known, non colliding indices//
    dnadb.insert(DNA("ACGTT", RndLocation.getRandNum(), true));
    //should be index 70//
    dnadb.insert(DNA("CGTAT", RndLocation.getRandNum(), true));
    //should be 49
    dnadb.insert(DNA("TTAGA", RndLocation.getRandNum(), true));
    //should be 42//
    if(dnadb.m_currentTable[70]->m_sequence == "ACGTT" && dnadb.m_currentTable[49]->m_sequence == "CGTAT" && dnadb.m_currentTable[42]->m_sequence == "TTAGA" && dnadb.m_currentSize == 13){
        return true;
    }
    return false;
}
//error case for getdna//
bool Tester::getDNAErrorCase(){
    Random RndLocation(MINLOCID,MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode, LINEAR);
    for (int i=0;i<10;i++){
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum(), true);
        dnadb.insert(dataObj);
    }
    //doesnt exist in database, but searches for it//
    DNA found = dnadb.getDNA("APPLE", 123456);
    //true if emty object is returned//
    if(found.m_sequence == "" && found.m_location == 0){
        return true;
    }
    return false;
}
//standard case for getdna//
bool Tester::getDNAstandardcase(){
    Random RndLocation(MINLOCID,MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode, LINEAR);
    for (int i=0;i<40;i++){
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum(), true);
        dnadb.insert(dataObj);
    }
    //inserts a known object and searches for it//
    dnadb.insert(DNA("ACGTT", 123456, true));
    DNA found = dnadb.getDNA("ACGTT", 123456);
    if(found.m_sequence == "ACGTT" && found.m_location == 123456){
        return true;
    }
    return false;
}
//collision case for getdna//
bool Tester::getDNAcollisioncase(){
    Random RndLocation(MINLOCID,MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode, LINEAR);
    for (int i=0;i<40;i++){
        // generating random data
        DNA dataObj = DNA("ACGTT", RndLocation.getRandNum(), true);
        dnadb.insert(dataObj);
    }
    //inserts multiple objects with the same sequence to cause collisions//
    dnadb.insert(DNA("ACGTT", 123456, true));
    //searches for one specifically and verifies the location//
    DNA found = dnadb.getDNA("ACGTT", 123456);
    if(found.m_sequence == "ACGTT" && found.m_location == 123456){
        return true;
    }
    return false;
}

//standard case for remove function//
bool Tester::removeTestStandardCase(){
    Random RndLocation(MINLOCID,MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode, QUADRATIC);
    for (int i=0;i<40;i++){
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum(), true);
        dnadb.insert(dataObj);
    }
    //inserts objects into table, removes known object//
    dnadb.insert(DNA("ACGTT", 123456, true));
    DNA node("ACGTT", 123456, true);
    if(dnadb.remove(node)){
        return true;
    }
    return false;
}
//tests if rehash is triggered if load factor hits capacity//
bool Tester::rehashTest1(){
    vector<DNA> datalist;
    Random RndLocation(MINLOCID,MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode, LINEAR);
    for (int i=0;i<51;i++){
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum(), true);
        if(dnadb.insert(dataObj)){
            datalist.push_back(dataObj);
        }
    }
    //fails if the old table isnt created//
    if(dnadb.m_oldTable == nullptr){
        return false;
    }
    bool result = true;
    //since old table is created, checks if all data is still present in both tables//
    for (vector<DNA>::iterator it = datalist.begin(); it != datalist.end(); it++){
        DNA anObj = dnadb.getDNA((*it).getSequence(), (*it).getLocId());
        bool foundIt = (*it == anObj);
        result = result && foundIt;
        if (!foundIt){
            cout << "Data point " << 
            (*it).getSequence() << "(" <<
            (*it).getLocId() << ")" <<
            " is missing!" << endl;
        }
    }
    if(result){
        return true;
    }
    
}
//tests if rehash will complete//
bool Tester::rehashTest2(){
    vector<DNA> datalist;
    Random RndLocation(MINLOCID,MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode, LINEAR);
    for (int i=0;i<140;i++){
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum(), true);
        if(dnadb.insert(dataObj)){
            datalist.push_back(dataObj);
        }
    }
    bool result = true;
    //inserts enough objects to complete rehashing//
    for (vector<DNA>::iterator it = datalist.begin(); it != datalist.end(); it++){
        DNA anObj = dnadb.getDNA((*it).getSequence(), (*it).getLocId());
        bool foundIt = (*it == anObj);
        result = result && foundIt;
        if (!foundIt){
            cout << "Data point " << 
            (*it).getSequence() << "(" <<
            (*it).getLocId() << ")" <<
            " is missing!" << endl;
        }
    }
    //if all data points are present and the old table gets deleted, passes//
    if(result && dnadb.m_oldTable == nullptr){
        return true;
    }
    return false;

    
}
//tests if enough removals trigger a rehash//
bool Tester::rehashTest3() {
    Random RndLocation(MINLOCID, MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode, DOUBLEHASH);
    int originalSize = 0;
    for (int i = 0; i < 45; i++) {
        int location = RndLocation.getRandNum();
        DNA dataObj(sequencer(5, i), location, true);
        if(dnadb.insert(dataObj)){
            originalSize++;
        }
    }
    //inserts nodes into tabl;e//
    int removedCount = 0;
    //removes enough nodes to start a rehash//
    for(int i = 100; i > 49; i--){
        if(dnadb.m_currentTable[i] != nullptr && dnadb.m_currentTable[i] -> m_used){
            if(dnadb.remove(*dnadb.m_currentTable[i])){
                removedCount++;
            }
        }
    }
    //if the current size of the table plus the number of removed nodes equals the original size, returns true//
    if(dnadb.m_currentSize + removedCount == originalSize){
        return true;
    }

    return false;
}
//tests if rehash will complete if enough are removed//
bool Tester::rehashTest4() {
    Random RndLocation(MINLOCID, MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode, DOUBLEHASH);
    int originalSize = 0;
    for (int i = 0; i < 45; i++) {
        int location = RndLocation.getRandNum();
        DNA dataObj(sequencer(5, i), location, true);
        if(dnadb.insert(dataObj)){
            originalSize++;
        }
    }
    //same logic, inserts, then removes enough to trigger a rehash//
    int removedCount = 0;
    for(int i = 100; i > 0; i--){
        if(dnadb.m_currentTable[i] != nullptr && dnadb.m_currentTable[i] -> m_used){
            if(dnadb.remove(*dnadb.m_currentTable[i])){
                removedCount++;
            }
        }
    }
    //if all data is accounted for, returns true//
    if(dnadb.m_currentSize + removedCount == originalSize){
        return true;
    }

    return false;
}

//tests insertion function//
bool Tester::quadraticInsertionTest(){
    Random RndLocation(MINLOCID,MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode, QUADRATIC);
    for (int i=0;i<10;i++){
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum(), true);
        dnadb.insert(dataObj);
    }
    //adding data points with known, colliding indices to force a quadratic probe//
    dnadb.insert(DNA("ACGTT", 123456, true));
    //should be index 70//
    dnadb.insert(DNA("ACGTT", 654321, true));
    //should be 71//
    dnadb.insert(DNA("ACGTT", 333333, true));
    //should be 75//
    if(dnadb.m_currentTable[70]->m_sequence == "ACGTT" && 
    dnadb.m_currentTable[70] -> m_location == 123456 &&
    dnadb.m_currentTable[71]->m_sequence == "ACGTT" && 
    dnadb.m_currentTable[71] -> m_location == 654321 &&
    dnadb.m_currentTable[75]->m_sequence == "ACGTT" && 
    dnadb.m_currentTable[75] -> m_location == 333333 &&
    dnadb.m_currentSize == 13){
        return true;
    }
    return false;
}

//tests insertion function//
bool Tester::doublehashInsertionTest(){
    Random RndLocation(MINLOCID,MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode, DOUBLEHASH);
    for (int i=0;i<10;i++){
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum(), true);
        dnadb.insert(dataObj);
    }
    //adding data points with known, colliding indices to force a quadratic probe//
    dnadb.insert(DNA("ACGTT", 123456, true));
    //should be index 70//
    dnadb.insert(DNA("ACGTT", 654321, true));
    //should be 74//
    dnadb.insert(DNA("ACGTT", 333333, true));
    //should be 78//
    if(dnadb.m_currentTable[70]->m_sequence == "ACGTT" && 
    dnadb.m_currentTable[70] -> m_location == 123456 &&
    dnadb.m_currentTable[74]->m_sequence == "ACGTT" && 
    dnadb.m_currentTable[74] -> m_location == 654321 &&
    dnadb.m_currentTable[78]->m_sequence == "ACGTT" && 
    dnadb.m_currentTable[78] -> m_location == 333333 &&
    dnadb.m_currentSize == 13){
        return true;
    }
    return false;
}
//remove function error case//
bool Tester::removeErrorCase(){
    Random RndLocation(MINLOCID,MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode, DOUBLEHASH);
    for (int i=0;i<49;i++){
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum(), true);
        dnadb.insert(dataObj);
    }
    //attempting to remove a string that doesn't exist in the table//
    DNA node("APPLE", 111111, true);
    if(!dnadb.remove(node)){
        return true;
    }
    return false;
}
//update location normal case//
bool Tester::updateLocNormalCase(){
    Random RndLocation(MINLOCID,MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode, DOUBLEHASH);
    for (int i=0;i<30;i++){
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum(), true);
        dnadb.insert(dataObj);
    }
    //should be at index 15//
    dnadb.insert(DNA("TAAAG", 123456,true));
    dnadb.updateLocId(DNA("TAAAG", 123456,true),654321);
    if(dnadb.m_currentTable[15]->m_location == 654321){
        return true;
    }
    return false;
}

//update location normal case//
bool Tester::updateLocErrorCase(){
    Random RndLocation(MINLOCID,MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode, DOUBLEHASH);
    for (int i=0;i<30;i++){
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum(), true);
        dnadb.insert(dataObj);
    }
    //tries to update a nonexistent index//
    if(!dnadb.updateLocId(DNA("APPLES", 123456,true),654321)){
        return true;
    }
    return false;
}

unsigned int hashCode(const string str) {
    unsigned int val = 0 ;
    const unsigned int thirtyThree = 33 ;  // magic number from textbook
    for ( int i = 0 ; i < str.length(); i++)
       val = val * thirtyThree + str[i] ;
    return val ;
 }
 string sequencer(int size, int seedNum){
     //this function returns a random DNA sequence
     // size param specifies the size of string
     string sequence = "";
     Random rndObject(0,3);
     rndObject.setSeed(seedNum);
     for (int i=0;i<size;i++){
         sequence = sequence + ALPHA[rndObject.getRandNum()];
     }
     return sequence;
 }
