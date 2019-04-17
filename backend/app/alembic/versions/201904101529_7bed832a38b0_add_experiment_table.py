"""Add experiment table

Revision ID: 7bed832a38b0
Revises: b9935fae5e7c
Create Date: 2019-04-10 15:30:11.616899

"""
import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import JSONB

# revision identifiers, used by Alembic.

revision = '7bed832a38b0'
down_revision = 'b9935fae5e7c'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table('experiment',
                    sa.Column('id', sa.Integer(), primary_key=True, index=True),
                    sa.Column('name', sa.String(), nullable=False, index=True, unique=True),
                    sa.Column('description', sa.String()),
                    sa.Column('root_directory', sa.String(4096)),
                    sa.Column('location', sa.String(4096)),
                    sa.Column('meta', JSONB()),
                    sa.Column('created_at', sa.DateTime(), nullable=False),
                    )


def downgrade():
    op.drop_table('experiment')
