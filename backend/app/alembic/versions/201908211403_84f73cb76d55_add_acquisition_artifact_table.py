"""Add acquisition_artifact table

Revision ID: 84f73cb76d55
Revises: 35b87e175854
Create Date: 2019-08-21 14:03:14.835314

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB

# revision identifiers, used by Alembic.
revision = '84f73cb76d55'
down_revision = '35b87e175854'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'acquisition_artifact',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('acquisition_id', sa.Integer(), sa.ForeignKey("acquisition.id", ondelete="CASCADE"), index=True),
        sa.Column('type', sa.String(64), nullable=False, index=True),
        sa.Column('description', sa.String()),
        sa.Column('meta', JSONB()),
        sa.Column('location', sa.String(4096)),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
    )


def downgrade():
    op.drop_table('acquisition_artifact')
