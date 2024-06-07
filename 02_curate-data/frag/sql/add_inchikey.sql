ALTER TABLE molecules RENAME TO old_table;

CREATE TABLE molecules
(
id integer primary key,
smiles text unique,
inchikey text unique,
natoms int,
elements blob,
tag text
);

INSERT INTO molecules SELECT id, smiles, smiles, natoms, elements, tag FROM old_table;

DROP TABLE old_table;
