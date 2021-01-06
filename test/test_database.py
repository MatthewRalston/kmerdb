import sys
import os
import unittest
import tempfile

sys.path.append("..")

from kmerdb import database

class Test_parsefile(unittest.TestCase):

    def setUp(self):

        self.k = 3
        self.database  = tempfile.NamedTemporaryFile(mode="w+", suffix=".sqlite3", delete=False)
        self.database.close()


    def test_str_filepath_only(self):
        """
        SqliteKdb throws a TypeError if it receives anything but a str
        """
        with self.assertRaises(TypeError):
            database.SqliteKdb(None, self.k)
        with self.assertRaises(TypeError):
            database.SqliteKdb(1, self.k)
        with self.assertRaises(TypeError):
            database.SqliteKdb(1.0, self.k)
        with self.assertRaises(TypeError):
            database.SqliteKdb([1], self.k)
        with self.assertRaises(TypeError):
            database.SqliteKdb({'hello': 'world'}, self.k)

    def test_filepath_exists(self):
        """
        SqliteKdb throws an OSError if it receives a filepath that doesnt exist
        """
        with self.assertRaises(OSError):
            db = database.SqliteKdb("helloworld.txt", self.k)
            db.conn.close()
            os.unlink("helloworld.txt")

    def test_k_is_int_only(self):
        """
        SqliteKdb throws a TypeError if it receives anything but an int
        """
        with self.assertRaises(TypeError):
            database.SqliteKdb(self.database, None)
        with self.assertRaises(TypeError):
            database.SqliteKdb(self.database, True)
        with self.assertRaises(TypeError):
            database.SqliteKdb(self.database, 'hello')
        with self.assertRaises(TypeError):
            database.SqliteKdb(self.database, 1.0)
        with self.assertRaises(TypeError):
            database.SqliteKdb(self.database, [1])
        with self.assertRaises(TypeError):
            database.SqliteKdb(self.database, {'hello': 'world'})


    def test_returns_database(self):
        """
        SqliteKdb returns a database object
        """

        testDb = database.SqliteKdb(self.database.name, self.k)
        self.assertTrue(isinstance(testDb, database.SqliteKdb))


if __name__ == '__main__':
    unittest.main()
