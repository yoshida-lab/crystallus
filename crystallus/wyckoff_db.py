# Copyright 2021 TsumiNa
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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
