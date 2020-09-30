"""Initialization

Revision ID: 6f205a3fb0b0
Revises:
Create Date: 2020-09-30 11:15:52.529210

"""
from alembic import op


# revision identifiers, used by Alembic.
revision = '6f205a3fb0b0'
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    op.execute('CREATE EXTENSION IF NOT EXISTS "uuid-ossp";')  # For uuid_generate_v4


def downgrade():
    op.execute('DROP EXTENSION IF EXISTS "uuid-ossp";')
