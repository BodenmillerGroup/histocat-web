"""Add experiment table

Revision ID: 7bed832a38b0
Revises: b9935fae5e7c
Create Date: 2019-04-10 15:30:11.616899

"""
import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import JSONB, ARRAY

# revision identifiers, used by Alembic.

revision = '7bed832a38b0'
down_revision = 'b9935fae5e7c'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'experiment',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('user_id', sa.Integer(), index=True),
        sa.Column('name', sa.String(), nullable=False, index=True, unique=True),
        sa.Column('description', sa.String()),
        sa.Column('location', sa.String(4096)),
        sa.Column('meta', JSONB()),
        sa.Column('tags', ARRAY(sa.String(64)), index=True),
        sa.Column('created_at', sa.DateTime(), nullable=False),
    )
    op.create_foreign_key(
        'fk_experiment_user',
        'experiment', 'user',
        ['user_id'], ['id'],
    )


def downgrade():
    op.drop_constraint('fk_experiment_user')
    op.drop_table('experiment')
