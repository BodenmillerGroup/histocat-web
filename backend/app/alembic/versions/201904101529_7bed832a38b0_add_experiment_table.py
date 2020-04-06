"""Add experiment table

Revision ID: 7bed832a38b0
Revises: fed22b8a9121
Create Date: 2019-04-10 15:30:11.616899

"""
import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import JSONB, ARRAY

# revision identifiers, used by Alembic.

revision = '7bed832a38b0'
down_revision = 'fed22b8a9121'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'experiment',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('user_id', sa.Integer(), sa.ForeignKey("user.id", ondelete="CASCADE"), index=True),
        sa.Column('name', sa.String(), nullable=False, index=True, unique=True),
        sa.Column('description', sa.String()),
        sa.Column('meta', JSONB()),
        sa.Column('tags', ARRAY(sa.String(64)), index=True),
        sa.Column('is_public', sa.Boolean(), nullable=False, default=False),
        sa.Column('location', sa.String(4096)),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
    )


def downgrade():
    op.drop_table('experiment')
