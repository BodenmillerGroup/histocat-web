"""Add channel table

Revision ID: 96df6178591e
Revises: b182e0878258
Create Date: 2019-04-10 15:32:17.364696

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = '96df6178591e'
down_revision = 'b182e0878258'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table('channel',
                    sa.Column('id', sa.Integer(), primary_key=True, index=True),
                    sa.Column('acquisition_id', sa.Integer(), index=True),
                    sa.Column('name', sa.String(), nullable=False, index=True),
                    sa.Column('description', sa.String()),
                    sa.Column('location', sa.String(4096)),
                    sa.Column('meta', sa.JSON()),
                    sa.Column('created_at', sa.DateTime(), nullable=False),
                    )
    op.create_foreign_key(
        'fk_channel_acquisition',
        'channel', 'acquisition',
        ['acquisition_id'], ['id'],
    )


def downgrade():
    op.drop_constraint('fk_channel_acquisition')
    op.drop_table('channel')
