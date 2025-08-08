Uses a hashmap implementation and incremental rehashing to simulate a database full of DNA object entries. Includes an insert, search, and delete function as well as options to change the desired method of probing for searching.
The hash table begins incremental rehashing when a certain threshold is met, allocating memory for a new table and transferring 25% of the old table's data points after each insert/delete operation.
Utilizes unit testing for each major function.
