"""Add user table

Revision ID: 219260649256
Revises: 6f205a3fb0b0
Create Date: 2020-09-30 11:17:53.304009

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = '219260649256'
down_revision = '6f205a3fb0b0'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'user',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('name', sa.String(), index=True),
        sa.Column('email', sa.String(), unique=True, index=True, nullable=False),
        sa.Column('password', sa.String()),
        sa.Column('is_active', sa.Boolean(), nullable=False),
        sa.Column('is_admin', sa.Boolean(), nullable=False),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
        sa.Column('updated_at', sa.DateTime(), default=sa.sql.func.now(), onupdate=sa.sql.func.now(), nullable=False),
    )


def downgrade():
    op.drop_table('user')
