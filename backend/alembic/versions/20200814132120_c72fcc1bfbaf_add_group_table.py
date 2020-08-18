"""Add group table

Revision ID: c72fcc1bfbaf
Revises: a4d621b5f66c
Create Date: 2020-08-14 13:21:20.332666

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = 'c72fcc1bfbaf'
down_revision = 'a4d621b5f66c'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'group',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('name', sa.String(), index=True, unique=True, nullable=False),
        sa.Column('description', sa.String()),
        sa.Column('url', sa.String()),
        sa.Column('is_open', sa.Boolean(), default=False, nullable=False),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
    )


def downgrade():
    op.drop_table('group')
