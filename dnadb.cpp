// CMSC 341 - Spring 25 - Project 4
#include "dnadb.h"
//constructor//
DnaDb::DnaDb(int size, hash_fn hash, prob_t probing = DEFPOLCY){
    //sets cap based on given size and parameters//
    if((isPrime(size)) && size < MAXPRIME && size > MINPRIME){
        m_currentCap = size;
    }else if(size <= MINPRIME){
        m_currentCap = MINPRIME;
    }else if(size >= MAXPRIME){
        m_currentCap = MAXPRIME;
    }else if(!isPrime(size)){
        m_currentCap = findNextPrime(size);
    }
    //initializes variables and current table//
    m_currentSize = 0;
    m_hash = hash;
    m_currProbing = probing;
    m_newPolicy = DEFPOLCY;
    m_currNumDeleted = 0;
    m_transferIndex = 0;
    m_currentTable = new DNA*[m_currentCap];
    for(int i = 0; i < m_currentCap;i++){
        m_currentTable[i] = nullptr;
    }
    //old table kept nullptr until needed//
    m_oldTable = nullptr;
    m_oldCap = 0;
    m_oldNumDeleted = 0;
    m_oldSize = 0;
    m_oldProbing = probing;
}
//destructor//
DnaDb::~DnaDb(){
    //calls clear functions if needed//
    if(m_currentTable != nullptr){
        clear(m_currentTable,m_currentCap);
        m_currentTable = nullptr;
    }
    if(m_oldTable != nullptr){
        clear(m_oldTable,m_oldCap);
        m_oldTable = nullptr;
    }
    m_currentSize = 0;
    m_transferIndex = 0;
    m_currentCap = 0;
    m_currNumDeleted = 0;
    m_oldCap = 0;
    m_oldNumDeleted = 0;
    m_oldSize = 0;
}
//clear helper function//
void DnaDb::clear(DNA** table,int size){
    if(table != nullptr){
        for(int i = 0; i < size; i++){
            if(table[i] != nullptr){
                delete table[i];
                table[i] = nullptr;
            }
        }
        delete[] table;
    }
}
//sets probing policy//
void DnaDb::changeProbPolicy(prob_t policy){
    m_newPolicy = policy;
}
//insert function//
bool DnaDb::insert(DNA dna){
    //calculates initial index//
    int index = m_hash(dna.m_sequence) % m_currentCap;
    bool inserted = false;
    int insertindex = 0;
    bool already_exists = false;
    //if the sequence and location are already existing in the table, returns//
    for(int i = 0; i < m_currentCap;i++){
        if(m_currentTable[i] != nullptr && m_currentTable[i] -> m_sequence == dna.m_sequence && m_currentTable[i] -> m_location == dna.m_location){
            already_exists = true;
        }
    }
    if(already_exists){
        return false;
    }
    //if indx is currently unoccupied, creates new object at that index//
    if(m_currentTable[index] != nullptr && !m_currentTable[index] -> m_used){
        m_currentTable[index] = new DNA(dna.m_sequence,dna.m_location, true);
        inserted = true; 
    }else if(m_currProbing == LINEAR){
        //checks for duplicates with linear probing//
        if(linearsearch(dna.m_sequence,dna.m_location,index,m_currentTable,m_currentCap) != -1){
            return false;
        }
        //returns index using linear probing//
        insertindex = linearinsert(index,m_currentTable,m_currentCap);
        inserted = true;
    }else if(m_currProbing == QUADRATIC){
        if(quadraticsearch(dna.m_sequence,dna.m_location, index,m_currentTable,m_currentCap) != -1){
            return false;
        }
        //returns index with quadratic probing//
        insertindex = quadraticinsert(index,m_currentTable,m_currentCap);
        inserted = true;
    }else if(m_currProbing == DOUBLEHASH){
        if(doublesearch(dna.m_sequence, dna.m_location, index,m_currentTable,m_currentCap) != -1){
            return false;
        }
        //returns index using double hash probing//
        insertindex = doubleinsert(dna.m_sequence,index,m_currentTable,m_currentCap);
        inserted = true;
    }
    //if no index is found, returns false, otherwise, inserts at the found index and increments size//
    if(insertindex == -1){
        return false;
    }else{
        m_currentTable[insertindex] = new DNA(dna.m_sequence,dna.m_location, true);
        m_currentSize++;
    }
    //checks load factor and if needed, sets up rehashing//
    if((lambda() > 0.5)){
        rehashsetup();
    }
    //if oldtable already exists, begins rehashing//
    if(m_oldTable != nullptr){
        rehash();
    }
    return inserted;
}
//rehash setup helper//
void DnaDb::rehashsetup(){
    //guards against premature table swapping//
    if(m_oldTable != nullptr){
        return;
    }
    //gets new cap//
    int newsize = findNextPrime(m_currentSize * 4);
    int newcap;
    if (isPrime(newsize) && newsize < MAXPRIME && newsize > MINPRIME) {
        newcap = newsize;
    }else if (newsize < MINPRIME){
        newcap = MINPRIME;
    }else if (newsize > MAXPRIME){
        newcap = MAXPRIME;
    }else{
        newcap = findNextPrime(newsize);
    }
    //allocates and initializes new table and variables and then calls rehash//
    DNA** newtable = new DNA*[newcap];
    for(int i = 0; i < newcap; i++){
        newtable[i] = nullptr;
    }
    m_oldTable = m_currentTable;
    m_oldCap = m_currentCap;
    m_oldSize = m_currentSize;
    m_oldProbing = m_currProbing;
    m_oldNumDeleted = m_currNumDeleted;
    m_currentTable = newtable;
    m_currentCap = newcap;
    m_currentSize = 0;
    m_currNumDeleted = 0;
    if(m_newPolicy != DEFPOLCY){
        m_currProbing = m_newPolicy;
    }
    rehash();

}

//rehash helper//
void DnaDb::rehash(){
    //exists early if the old table doesnt exist or the transfer index is greater than or equal to the old cap//
    if(m_oldTable == nullptr || m_transferIndex >= m_oldCap){
        return;
    }
    //sets the end index for rehashing//
    int endindex = 0;
    if(m_transferIndex + (m_oldCap/4) > m_oldCap){
        endindex = m_oldCap;
    }else{
        endindex = m_transferIndex + (m_oldCap/4);
    }
    //scans 1/4th of the table at a time, rehashing active elements to the new table and deleting afterward//
    for(int i = m_transferIndex; i < endindex; i++){
        if(m_oldTable[i] != nullptr && m_oldTable[i]-> m_used){
            rehashinsert(m_oldTable[i] -> m_sequence, m_oldTable[i] -> m_location, m_oldTable[i] -> m_used);
            delete m_oldTable[i];
            m_oldTable[i] = nullptr;
            m_oldNumDeleted++;
        }else if(m_oldTable[i] != nullptr && !m_oldTable[i]-> m_used){
            delete m_oldTable[i];
            m_oldTable[i] = nullptr;
        }
    }
    //increments transfer index and checks if rehash is finished//
    m_transferIndex = endindex;
    if(m_transferIndex >= m_oldCap){
        clear(m_oldTable,m_oldCap);
        if(m_oldTable != nullptr){
            m_oldTable = nullptr;
        }
        m_oldSize = 0;
        m_transferIndex = 0;
    }
}
//rehash insert helper//
void DnaDb::rehashinsert(string sequence, int location, bool used){
    //gets the initial index//
    int index = m_hash(sequence) % m_currentCap;
    int insertindex = 0;
    //if the initial index is unoccupied or inactive, creates a new object//
    if(m_currentTable[index] != nullptr && !m_currentTable[index] -> m_used){
        delete m_currentTable[index];
        m_currentTable[index] = new DNA(sequence, location, true);
        return; 
    }else if(m_oldProbing == LINEAR){
        //uses linear,quadratic, and double hash helpers to place new data poitts
        if(linearsearch(sequence,location,index,m_currentTable,m_currentCap) != -1){
            return;
        }
        insertindex = linearinsert(index,m_currentTable,m_currentCap);
    }else if(m_oldProbing == QUADRATIC){
        if(quadraticsearch(sequence,location,index,m_currentTable,m_currentCap) != -1){
            return;
        }
        insertindex = quadraticinsert(index,m_currentTable,m_currentCap);
    }else if(m_oldProbing == DOUBLEHASH){
        if(doublesearch(sequence,location,index,m_currentTable,m_currentCap) != -1){
            return;
        }
        insertindex = doubleinsert(sequence,index,m_currentTable,m_currentCap);
    }
    //if a legitimate index is found, creates a new object//
    if(insertindex != -1){
        m_currentTable[insertindex] = new DNA(sequence,location,true);
        m_currentSize++;
        return;
    }
}

//linear insert//
int DnaDb::linearinsert(int index, DNA** table, int size){
    int currindex = index;
    int i = 0;
    //uses linear probing to find an open spot//
    while(i < size){
        if(table[currindex] == nullptr || !table[currindex] -> m_used){
            return currindex;
        }
        currindex = (currindex +1) % size;
        i++;
    }
    return -1;
}
//quadratic insert//
int DnaDb::quadraticinsert(int index, DNA** table, int size){
    int currindex = index;
    int i = 0;
    //uses quadratic probing formula to find an open spot//
    while(i < size){
        currindex = (currindex + (i*i)) % size;
        if(table[currindex] == nullptr || !table[currindex] -> m_used){
            return currindex;
        }
        i++;
    }
    return -1;

}
//double hash insert//
int DnaDb::doubleinsert(string sequence, int index, DNA** table, int size){
    int currindex = index;
    int i = 0;
    int hashedsequence = m_hash(sequence);
    //uses the given double hash formula to hash object//
    while(i < size){
        currindex = ((hashedsequence % size) + i * (11-(hashedsequence % 11))) % size;
        if(table[currindex] == nullptr || !table[currindex] -> m_used){
            return currindex;
        }
        i++;
    }
    return -1;
}

//linear search//
int DnaDb::linearsearch(string sequence, int location, int index, DNA** table, int size)const{
    int currindex = index;
    int i = 0;
    while(i < size){
        //uses linear probing to find where an pbject should be//
        if(table[currindex] != nullptr && table[currindex] -> m_used){
            if(table[currindex] -> m_sequence == sequence && table[currindex] -> m_location == location){
                return currindex;
            }
        }
        currindex = (currindex +1) % size;
        i++;
    }
    return -1;
}
//quadratic search//
int DnaDb::quadraticsearch(string sequence, int location, int index, DNA** table, int size)const{
    int currindex = index;
    int i = 0;
    while(i < size){
        currindex = (currindex + (i*i)) % size;
        //same logic as quadratoc insert//
        if(table[currindex] != nullptr && table[currindex] -> m_used){
            if(table[currindex] -> m_sequence == sequence && table[currindex] -> m_location == location){
                return currindex;
            }
        }
        i++;
    }
    return -1;

}
//double hash search//
int DnaDb::doublesearch(string sequence, int location, int index, DNA** table, int size)const{
    int currindex = index;
    int i = 0;
    int hashedsequence = m_hash(sequence);
    while(i < size){
        //same logic as double insert//
        currindex = ((hashedsequence % size) + i * (11-(hashedsequence % 11))) % size;
        if(table[currindex] != nullptr && table[currindex] -> m_used){
            if(table[currindex] -> m_sequence == sequence && table[currindex] -> m_location == location){
                return currindex;
            }
        }
        i++;
    }
    return -1;
}

//remove function//
bool DnaDb::remove(DNA dna){
    //gets initial hash//
    int index = m_hash(dna.m_sequence) % m_currentCap;
    int oldindex = 1;
    if(m_oldTable != nullptr && m_oldCap != 0){
        oldindex = m_hash(dna.m_sequence) % m_oldCap;
    }
    int removedindex = -1;
    bool old = false;
    //if the initial index is the desired node, tags it as deleted//
    if(m_currentTable[index] != nullptr && 
       m_currentTable[index]->m_sequence == dna.m_sequence &&
       m_currentTable[index]->m_location == dna.m_location &&
       m_currentTable[index]->m_used){
        m_currentTable[index]->m_used = false;
        m_currNumDeleted++;
        m_currentSize--;
        return true;
    }
    //if not found immediately, uses probing helpers to find index//
    if(m_currProbing == LINEAR){
        removedindex = linearsearch(dna.m_sequence,dna.m_location,index,m_currentTable,m_currentCap); 
        //if it wasn't found in the current table, checks the old table//
        if(removedindex == -1 && m_oldTable != nullptr){
            removedindex = linearsearch(dna.m_sequence,dna.m_location,oldindex,m_oldTable,m_oldCap);
            if(removedindex != -1){
                old = true;
            }
        }
    }else if(m_currProbing == QUADRATIC){
        removedindex = quadraticsearch(dna.m_sequence,dna.m_location,index,m_currentTable,m_currentCap); 
        if(removedindex == -1 && m_oldTable != nullptr){
            removedindex = quadraticsearch(dna.m_sequence,dna.m_location,oldindex,m_oldTable,m_oldCap);
            if(removedindex != -1){
                old = true;
            }
        }
    }else if(m_currProbing == DOUBLEHASH){
        removedindex = doublesearch(dna.m_sequence,dna.m_location,index,m_currentTable,m_currentCap); 
        if(removedindex == -1 && m_oldTable != nullptr){
            removedindex = doublesearch(dna.m_sequence,dna.m_location,oldindex,m_oldTable,m_oldCap);
            if(removedindex != -1){
                old = true;
            }
        }
    }
    //if not found, returns false//
    if(removedindex == -1){
        return false;
    }
    //if old flag is true, deletes from old table, if not, deletes from cureent table//
    if(old){
        m_oldTable[removedindex]->m_used = false;
        m_oldNumDeleted++;
        m_oldSize--;
    }else{
        m_currentTable[removedindex]->m_used = false;
        m_currNumDeleted++;
        m_currentSize--;
    }
    //checks deleted ratio if rehash is needed//
    if((deletedRatio() > 0.8)){
        rehashsetup();
    }
    //if rehashing is happening, calls rehash//
    if(m_oldTable != nullptr){
        rehash();
    }

    return true;
}

//getdna function//
const DNA DnaDb::getDNA(string sequence, int location) const{
    //gets initial hash//
    int index = m_hash(sequence) % m_currentCap;
    int oldindex = 0;
    //if the old table exists, gets the hash from the old table size//
    if(m_oldTable != nullptr && m_oldSize != 0){
        oldindex = m_hash(sequence) % m_oldCap;
    }
    bool old = false;
    //if the object is found immediately, returns an object//
    if(m_currentTable[index] != nullptr && m_currentTable[index] != nullptr && m_currentTable[index] -> m_sequence == sequence && m_currentTable[index] -> m_location == location && m_currentTable[index]-> m_used){
        return DNA(m_currentTable[index] -> m_sequence, m_currentTable[index] -> m_location, m_currentTable[index] -> m_used);
    }else{
        //if not, uses probing methods to find index//
        int newindex = 0;
        if(m_currProbing == LINEAR){
            newindex = linearsearch(sequence,location,index,m_currentTable,m_currentCap);
            if(newindex == -1){
                //checks old if not found in new//
                if(m_oldTable != nullptr){
                    newindex = linearsearch(sequence,location,oldindex,m_oldTable,m_oldCap);
                    if(newindex != -1){
                        old = true;
                    }
                }
            }
        }else if(m_currProbing == QUADRATIC){
            newindex = quadraticsearch(sequence,location,index,m_currentTable,m_currentCap);
            if(newindex == -1){
                if(m_oldTable != nullptr){
                    newindex = quadraticsearch(sequence,location,oldindex,m_oldTable,m_oldCap);
                    if(newindex != -1){
                        old = true;
                    }
                }
            }
        }else{
            newindex = doublesearch(sequence,location,index,m_currentTable,m_currentCap);
            if(newindex == -1){
                if(m_oldTable != nullptr){
                    newindex = doublesearch(sequence,location,oldindex,m_oldTable,m_oldCap);
                    if(newindex != -1){
                        old = true;
                    }
                }
            }
        }
        //if not found, returns an empty object//
        if(newindex == -1){
            return DNA();
        }else{
            //if in old, returns from old, else returns from new//
            if(old){
                return DNA(m_oldTable[newindex] -> m_sequence, m_oldTable[newindex] -> m_location, m_oldTable[newindex] -> m_used);
            }else{
                return DNA(m_currentTable[newindex] -> m_sequence, m_currentTable[newindex] -> m_location, m_currentTable[newindex] -> m_used);
            }
        }
    }
}

bool DnaDb::updateLocId(DNA dna, int location){
    //first checks the current table for the sequence, changes location and returns if found//
    for(int i = 0; i < m_currentCap; i++){
        if(m_currentTable[i] != nullptr && m_currentTable[i]->m_sequence == dna.m_sequence){
            m_currentTable[i] -> m_location = location;
            return true;
        }
    }
    //if not in curret table, checks the old table//
    if(m_oldTable != nullptr){
        for(int i = 0; i < m_oldCap; i++){
            if(m_oldTable[i] != nullptr && m_oldTable[i]->m_sequence == dna.m_sequence){
                m_oldTable[i] -> m_location = location;
                return true;
            }
        }
    }    
    //if not found, returns false//
    return false;
}
//calculates load factor//
float DnaDb::lambda() const {
    float loadfactor =  float(m_currentSize) / m_currentCap;
    return loadfactor;
}
//calculates deleted ratio//
float DnaDb::deletedRatio() const {
    float ratio = float(m_currNumDeleted)/m_currentSize;
    return ratio;
}

void DnaDb::dump() const {
    cout << "Dump for the current table: " << endl;
    if (m_currentTable != nullptr){
        for (int i = 0; i < m_currentCap; i++) {
            cout << "[" << i << "] : " << m_currentTable[i] << endl;
        }
    }
    cout << "Dump for the old table: " << endl;
    if (m_oldTable != nullptr){
        for (int i = 0; i < m_oldCap; i++) {
            cout << "[" << i << "] : " << m_oldTable[i] << endl;
        }
    }
}

bool DnaDb::isPrime(int number){
    bool result = true;
    for (int i = 2; i <= number / 2; ++i) {
        if (number % i == 0) {
            result = false;
            break;
        }
    }
    return result;
}

int DnaDb::findNextPrime(int current){
    //we always stay within the range [MINPRIME-MAXPRIME]
    //the smallest prime starts at MINPRIME
    if (current < MINPRIME) current = MINPRIME-1;
    for (int i=current; i<MAXPRIME; i++) { 
        for (int j=2; j*j<=i; j++) {
            if (i % j == 0) 
                break;
            else if (j+1 > sqrt(i) && i != current) {
                return i;
            }
        }
    }
    //if a user tries to go over MAXPRIME
    return MAXPRIME;
}
