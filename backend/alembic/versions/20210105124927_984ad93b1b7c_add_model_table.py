"""Add model table

Revision ID: 984ad93b1b7c
Revises: bcd5ecf42c3e
Create Date: 2021-01-05 12:49:27.310783

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB


# revision identifiers, used by Alembic.
revision = '984ad93b1b7c'
down_revision = 'bcd5ecf42c3e'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'model',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('name', sa.String(), nullable=False, index=True, unique=True),
        sa.Column('application', sa.String(), nullable=False),
        sa.Column('description', sa.String()),
        sa.Column('location', sa.String(4096)),
        sa.Column('meta', JSONB()),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
    )


def downgrade():
    op.drop_table('model')
