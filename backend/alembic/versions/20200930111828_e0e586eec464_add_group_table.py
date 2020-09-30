"""Add group table

Revision ID: e0e586eec464
Revises: 219260649256
Create Date: 2020-09-30 11:18:28.563462

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import ARRAY


# revision identifiers, used by Alembic.
revision = 'e0e586eec464'
down_revision = '219260649256'
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
        sa.Column('tags', ARRAY(sa.String(64)), index=True),
        sa.Column('location', sa.String(4096)),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
    )


def downgrade():
    op.drop_table('group')
