"""Add channel table

Revision ID: 96df6178591e
Revises: b182e0878258
Create Date: 2019-04-10 15:56:17.364696

"""
import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import JSONB

# revision identifiers, used by Alembic.
revision = '96df6178591e'
down_revision = 'b182e0878258'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'channel',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('acquisition_id', sa.Integer(), sa.ForeignKey("acquisition.id", ondelete="CASCADE"), index=True),
        sa.Column('metaname', sa.String(4096), index=True),
        sa.Column('original_id', sa.Integer(), index=True),
        sa.Column('metal', sa.String(), index=True),
        sa.Column('label', sa.String(), index=True),
        sa.Column('mass', sa.Integer()),
        sa.Column('max_intensity', sa.Integer()),
        sa.Column('min_intensity', sa.Integer()),
        sa.Column('meta', JSONB()),
        sa.Column('location', sa.String(4096)),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
    )


def downgrade():
    op.drop_table('channel')
