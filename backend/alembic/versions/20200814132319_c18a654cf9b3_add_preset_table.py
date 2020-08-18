"""Add preset table

Revision ID: c18a654cf9b3
Revises: 601250049059
Create Date: 2020-08-14 13:23:19.279740

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB


# revision identifiers, used by Alembic.
revision = 'c18a654cf9b3'
down_revision = '601250049059'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'preset',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('experiment_id', sa.Integer(), sa.ForeignKey("experiment.id", ondelete="CASCADE"), nullable=False, index=True),
        sa.Column('name', sa.String()),
        sa.Column('description', sa.String()),
        sa.Column('data', JSONB()),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
    )


def downgrade():
    op.drop_table('preset')
