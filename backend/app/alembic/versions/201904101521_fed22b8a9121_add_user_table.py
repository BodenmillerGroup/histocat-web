"""Add user table

Revision ID: fed22b8a9121
Revises: b9935fae5e7c
Create Date: 2019-04-10 15:29:17.676102

"""
import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.

revision = 'fed22b8a9121'
down_revision = 'b9935fae5e7c'
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
