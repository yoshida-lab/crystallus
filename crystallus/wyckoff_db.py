from peewee import Model, SqliteDatabase, CharField, ForeignKeyField, IntegerField, TextField, BooleanField
import os

db_path = os.path.dirname(os.path.abspath(__file__)) + '/spacegroup.db'
db = SqliteDatabase(db_path)


class BaseModel(Model):

    class Meta:
        database = db


class SpaceGroup(BaseModel):
    full_symbol = CharField(20, unique=True)
    symbol = CharField(20)
    symmetry = CharField(20)
    point_group = CharField(20)
    patterson_symmetry = CharField(20)
    spacegroup_num = IntegerField()
    hm_num = IntegerField()
    wy_sets = IntegerField()


class Wyckoff(BaseModel):
    space_group = ForeignKeyField(SpaceGroup, backref='wyckoffs')
    multiplicity = IntegerField()
    site_symmetry = CharField(20)
    letter = CharField(1)
    positions = TextField()
    reuse = BooleanField()
