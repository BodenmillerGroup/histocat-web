"""Initialization

Revision ID: a4d621b5f66c
Revises:
Create Date: 2020-08-14 13:21:00.329595

"""
from alembic import op


# revision identifiers, used by Alembic.
revision = 'a4d621b5f66c'
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    op.execute('CREATE EXTENSION IF NOT EXISTS "uuid-ossp";')  # For uuid_generate_v4


def downgrade():
    op.execute('DROP EXTENSION IF EXISTS "uuid-ossp";')

