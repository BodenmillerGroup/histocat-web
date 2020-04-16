"""Initialization

Revision ID: b9935fae5e7c
Revises:
Create Date: 2019-04-10 15:20:17.676102

"""
from alembic import op

# revision identifiers, used by Alembic.

revision = 'b9935fae5e7c'
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    op.execute('CREATE EXTENSION IF NOT EXISTS "uuid-ossp";')  # For uuid_generate_v4


def downgrade():
    op.execute('DROP EXTENSION IF EXISTS "uuid-ossp";')
