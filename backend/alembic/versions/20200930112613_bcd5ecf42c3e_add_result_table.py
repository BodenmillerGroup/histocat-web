"""Add result table

Revision ID: bcd5ecf42c3e
Revises: cf4962b43209
Create Date: 2020-09-30 11:26:13.381875

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB


# revision identifiers, used by Alembic.
revision = 'bcd5ecf42c3e'
down_revision = 'cf4962b43209'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'result',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('dataset_id', sa.Integer(), sa.ForeignKey("dataset.id", ondelete="CASCADE"), index=True, nullable=False),
        sa.Column('status', sa.String(64), index=True, nullable=False),
        sa.Column('name', sa.String()),
        sa.Column('description', sa.String()),
        sa.Column('pipeline', JSONB()),
        sa.Column('input', JSONB()),
        sa.Column('output', JSONB()),
        sa.Column('location', sa.String(4096)),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
    )


def downgrade():
    op.drop_table('result')
